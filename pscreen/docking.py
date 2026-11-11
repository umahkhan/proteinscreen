"""GNINA docking driver with a resumable SQLite job store.

One job = (ligand, pocket, mode). Finished jobs are skipped on re-run; failures are recorded and never abort
the run. Scores are stored as GNINA reports them and are NOT combined:
  minimizedAffinity  Vina-style score, kcal/mol (more negative = better)
  CNNscore           pose-quality probability, 0-1 (higher = more likely a good pose)
  CNNaffinity        predicted pK affinity (higher = tighter), CNN_VS = CNNscore x CNNaffinity
"""
from __future__ import annotations
import hashlib, json, os, sqlite3, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict, field
from pathlib import Path
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parent.parent
GNINA_SH = ROOT / "infra" / ("gnina.sh" if sys.platform == "darwin" else "gnina_linux.sh")   # macOS: Lima VM; Linux: native binary

MODE_NONCOVALENT = "noncovalent"
MODE_COVALENT = "covalent"

POCKET_SOURCES = ("user_box", "known_ligand", "predicted_conventional", "predicted_cryptic")


@dataclass
class Pocket:
    pocket_id: str
    center: tuple[float, float, float]
    size: tuple[float, float, float] = (22.0, 22.0, 22.0)
    source: str = "user_box"            # see POCKET_SOURCES
    receptor_path: str = ""             # conformation the pocket belongs to (cryptic pockets: the open/holo-like one)
    residues: list[str] = field(default_factory=list)   # e.g. ["A:ASP25", ...]
    metrics: dict = field(default_factory=dict)         # tool-specific scores, never merged
    @property
    def experimentally_confirmed(self) -> bool:
        return self.source == "known_ligand"


@dataclass
class DockConfig:
    exhaustiveness: int = 8
    num_modes: int = 9
    cpu: int = 4
    seed: int = 0
    cnn_scoring: str = "rescore"       # none|rescore|refinement|all  (rescore = CNN on final poses)
    autobox_add: float = 0.0
    timeout_s: int = 3600
    workers: int = 1                   # concurrent gnina processes (workers*cpu <= VM cores)
    def hash(self) -> str:
        return hashlib.sha1(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:10]


def gnina_version() -> str:
    out = subprocess.run([str(GNINA_SH), "--version"], capture_output=True, text=True, timeout=120)
    return (out.stdout or out.stderr).strip().splitlines()[0]


class Store:
    def __init__(self, path: str | Path):
        self.path = str(path)
        self.db = sqlite3.connect(self.path, check_same_thread=False)
        self.db.executescript("""
        CREATE TABLE IF NOT EXISTS jobs(ligand_key TEXT, compound_id TEXT, pocket_id TEXT, mode TEXT,
            status TEXT, error TEXT, seconds REAL, config_hash TEXT, finished_at TEXT,
            PRIMARY KEY(ligand_key, pocket_id, mode));
        CREATE TABLE IF NOT EXISTS poses(ligand_key TEXT, pocket_id TEXT, mode TEXT, pose_idx INTEGER,
            minimized_affinity REAL, cnn_score REAL, cnn_affinity REAL, cnn_vs REAL, cnn_affinity_variance REAL,
            pose_sdf TEXT, PRIMARY KEY(ligand_key, pocket_id, mode, pose_idx));
        CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
        """)
    def done(self, key, pocket, mode) -> bool:
        r = self.db.execute("SELECT status FROM jobs WHERE ligand_key=? AND pocket_id=? AND mode=?", (key, pocket, mode)).fetchone()
        return bool(r and r[0] == "done")
    def put_meta(self, k, v): self.db.execute("INSERT OR REPLACE INTO meta VALUES(?,?)", (k, json.dumps(v))); self.db.commit()
    def record(self, job, status, error, seconds, cfg_hash, poses):
        with self.db:
            self.db.execute("DELETE FROM poses WHERE ligand_key=? AND pocket_id=? AND mode=?", (job["ligand_key"], job["pocket_id"], job["mode"]))
            self.db.execute("INSERT OR REPLACE INTO jobs VALUES(?,?,?,?,?,?,?,?,datetime('now'))",
                (job["ligand_key"], job["compound_id"], job["pocket_id"], job["mode"], status, error, seconds, cfg_hash))
            self.db.executemany("INSERT INTO poses VALUES(?,?,?,?,?,?,?,?,?,?)", poses)
    def jobs(self) -> pd.DataFrame: return pd.read_sql("SELECT * FROM jobs", self.db)
    def poses(self) -> pd.DataFrame: return pd.read_sql("SELECT * FROM poses", self.db)


def _parse_poses(sdf: Path, job) -> list[tuple]:
    rows = []
    for i, m in enumerate(Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False)):
        if m is None: continue
        g = lambda k: float(m.GetProp(k)) if m.HasProp(k) else None
        rows.append((job["ligand_key"], job["pocket_id"], job["mode"], i, g("minimizedAffinity"), g("CNNscore"),
                     g("CNNaffinity"), g("CNN_VS"), g("CNNaffinity_variance"), str(sdf)))
    return rows


def _run_job(job, receptor: Path, pocket: Pocket, cfg: DockConfig, posedir: Path):
    out = posedir / f"{job['ligand_key']}__{pocket.pocket_id}__{job['mode']}.sdf"
    cmd = [str(GNINA_SH), "-r", str(receptor), "-l", job["ligand_sdf"], "-o", str(out),
           "--center_x", f"{pocket.center[0]:.3f}", "--center_y", f"{pocket.center[1]:.3f}", "--center_z", f"{pocket.center[2]:.3f}",
           "--size_x", f"{pocket.size[0]:.2f}", "--size_y", f"{pocket.size[1]:.2f}", "--size_z", f"{pocket.size[2]:.2f}",
           "--exhaustiveness", str(cfg.exhaustiveness), "--num_modes", str(cfg.num_modes), "--cpu", str(cfg.cpu),
           "--seed", str(cfg.seed), "--cnn_scoring", cfg.cnn_scoring] + ([] if os.environ.get("PSCREEN_GPU") == "1" else ["--no_gpu"])
    cov = job.get("covalent")
    if cov:   # GNINA covalent docking: ligand must already be in bound-state chemistry (see covalent.py)
        cmd += ["--covalent_rec_atom", cov["rec_atom"], "--covalent_lig_atom_pattern", cov["lig_smarts"],
                "--covalent_bond_order", str(cov.get("bond_order", 1))]
    t0 = time.time()
    env = {**os.environ, "GNINA_TIMEOUT": str(cfg.timeout_s)}
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=cfg.timeout_s + 60)
    except subprocess.TimeoutExpired:
        return job, "failed", f"timeout after {cfg.timeout_s}s", time.time() - t0, []
    dt = time.time() - t0
    if p.returncode != 0 or not out.exists():
        tail = (p.stderr or p.stdout).strip().splitlines()[-3:]
        reason = "gnina timed out" if p.returncode in (124, 137) else "gnina error: " + " | ".join(tail)
        return job, "failed", reason, dt, []
    poses = _parse_poses(out, job)
    if not poses:
        return job, "failed", "gnina produced no poses (ligand could not be placed in the search box)", dt, []
    return job, "done", None, dt, poses


def run_screen(prepared: pd.DataFrame, receptor: str | Path, pockets: list[Pocket], cfg: DockConfig,
               outdir: str | Path, progress=None, target_residue: str | None = None,
               pairs: set[tuple[str, str]] | None = None, mode: str = MODE_NONCOVALENT,
               covalent_specs: dict[tuple[str, str], dict] | None = None) -> Store:
    """Dock all status=='ok' ligands into all pockets (or only `pairs` of (ligand_key, pocket_id)). Resumable: call again with the same outdir to continue."""
    outdir = Path(outdir).resolve(); posedir = outdir / "poses"; posedir.mkdir(parents=True, exist_ok=True)
    store = Store(outdir / "results.sqlite")
    store.put_meta("config", asdict(cfg)); store.put_meta("pockets", [asdict(p) for p in pockets])
    store.put_meta("target_residue", target_residue or "")
    store.put_meta("gnina_version", gnina_version()); store.put_meta("receptor", str(receptor))
    if "ligand_sdf" not in prepared or not (prepared["status"] == "ok").any():
        raise ValueError("No dockable compounds: every record is invalid, duplicate or failed preparation (see prepared.csv for reasons).")
    ok = prepared[(prepared["status"] == "ok") & prepared["ligand_sdf"].notna()]
    todo = []
    for _, r in ok.iterrows():
        for pk in pockets:
            if pairs is not None and (r["ligand_key"], pk.pocket_id) not in pairs:
                continue
            if mode == MODE_COVALENT and (r["ligand_key"], pk.pocket_id) not in (covalent_specs or {}):
                continue
            if not store.done(r["ligand_key"], pk.pocket_id, mode):
                todo.append(({"ligand_key": r["ligand_key"], "compound_id": r["compound_id"], "ligand_sdf": str(Path(r["ligand_sdf"]).resolve()),
                              "pocket_id": pk.pocket_id, "mode": mode,
                              "covalent": (covalent_specs or {}).get((r["ligand_key"], pk.pocket_id)) if mode == MODE_COVALENT else None}, pk))
    total, n = len(todo), 0
    by_id = {p.pocket_id: p for p in pockets}
    with ThreadPoolExecutor(max_workers=cfg.workers) as ex:
        futs = [ex.submit(_run_job, j, Path(by_id[j["pocket_id"]].receptor_path or receptor).resolve(), pk, cfg, posedir) for j, pk in todo]
        for f in as_completed(futs):
            job, status, err, dt, poses = f.result()
            store.record(job, status, err, dt, cfg.hash(), poses)
            n += 1
            if progress: progress(n, total, job["compound_id"], status)
    return store
