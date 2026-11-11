"""End-to-end screening pipeline: prepare -> pockets -> GNINA docking -> covalent docking -> ADMET -> results/exports.

Run from the UI (background subprocess) or the CLI:   python -m pscreen.pipeline runs/<name>     (reads runs/<name>/run_config.json)
Every stage is idempotent and the docking stages are resumable (finished jobs are skipped), so re-running the same directory continues
a failed or interrupted run. Progress is written to progress.json; a reproducibility manifest (tool versions, input hashes, config) to
run_manifest.json. Failures in optional stages (pocket predictors, ADMET, covalent) are recorded and do not stop the run.
"""
from __future__ import annotations
import hashlib, json, platform, shutil, sys, time, traceback, zipfile
from dataclasses import dataclass, field, asdict
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from . import __version__
from .prep import prepare_library, summarize as prep_summary, embed_3d
from .docking import DockConfig, Pocket, run_screen, gnina_version, MODE_COVALENT
from .pockets import from_known_ligand, from_user_box
from .results import build_results, best_pose_per_compound
from .covalent import covalent_candidates
from .geometry import NUCLEOPHILE_ATOMS, load_structure, residue_atoms, parse_residue
from . import chem_props, hosting
from .warheads import get_warhead



@dataclass
class RunConfig:
    name: str
    receptor: str                      # path (PDB or mmCIF) as supplied
    library: str                       # path (.smi/.csv/.tsv/.sdf)
    target_residue: str = ""           # e.g. "A:CYS145" (primary; used for geometry columns)
    target_residues: list = field(default_factory=list)   # additional residues for the covalent lane (each is docked covalently)
    known_ligand: str = ""             # optional PDB/SDF/MOL2 defining the pocket (experimental reference)
    box_center: list = field(default_factory=list)   # optional [x,y,z]
    box_size: list = field(default_factory=lambda: [22.0, 22.0, 22.0])
    predict_pockets: bool = False
    predictors: list = field(default_factory=lambda: ["p2rank"])   # subset of p2rank|pocketminer|lacuna
    max_pockets_per_tool: int = 3
    pocketminer_threshold: float = 0.7
    box_padding: float = 6.0
    exhaustiveness: int = 8
    num_modes: int = 9
    cpu_per_job: int = 3
    workers: int = 2
    timeout_s: int = 1800
    seed: int = 0
    cnn_scoring: str = "rescore"
    covalent: bool = True
    covalent_max_center_A: float = 0.0   # 0 = dock covalently at every pocket (unbiased). >0 skips pockets whose centre is farther than this from the nucleophile (speed only; biases toward the target site)
    chemprop: bool = False
    chemprop_endpoints: list = field(default_factory=lambda: list(chem_props.DEFAULT_ENDPOINTS))
    experimental_csv: str = ""         # optional: compound_id + measured columns

    def dock_cfg(self) -> DockConfig:
        return DockConfig(self.exhaustiveness, self.num_modes, self.cpu_per_job, self.seed, self.cnn_scoring, 0.0, self.timeout_s, self.workers)

    def save(self, rundir: Path): (rundir / "run_config.json").write_text(json.dumps(asdict(self), indent=1))
    @staticmethod
    def load(rundir: Path) -> "RunConfig":
        d = json.loads((rundir / "run_config.json").read_text())
        return RunConfig(**{k: v for k, v in d.items() if k in RunConfig.__dataclass_fields__})

    def residues(self) -> list[str]:
        return list(dict.fromkeys([r for r in [self.target_residue, *self.target_residues] if r]))


STAGE_LABELS = {"inputs": "Validate inputs", "prepare": "Prepare compounds (RDKit, warheads)", "pockets": "Define binding sites", "docking": "GNINA docking (non-covalent)",
                "covalent_docking": "GNINA covalent docking", "properties": "ADMET-AI properties", "results": "Assemble results"}


def stage_label(key: str) -> str:
    if key.startswith("pocket:"): return {"p2rank": "P2Rank pockets", "pocketminer": "PocketMiner cryptic residues", "lacuna": "Lacuna ensemble pockets"}.get(key[7:], key[7:])
    return STAGE_LABELS.get(key, key)


def plan(cfg: "RunConfig") -> list[str]:
    """Ordered stage keys this configuration will run (the log panel shows them as pending before they start)."""
    st = ["inputs", "prepare", "pockets"]
    if cfg.predict_pockets or not (cfg.known_ligand or cfg.box_center): st += [f"pocket:{t}" for t in cfg.predictors]
    st.append("docking")
    if cfg.covalent and bool(cfg.residues()): st.append("covalent_docking")
    if cfg.chemprop: st.append("properties")
    return st + ["results"]


class Progress:
    """progress.json: per-stage state/timing (start, end epoch seconds), counters and messages, plus the planned stage order and run start/end."""
    def __init__(self, rundir: Path, stages: list[str] | None = None):
        self.p = rundir / "progress.json"; self.log = rundir / "run.log"
        self.state = {"state": "running", "started": time.time(), "plan": stages or [], "stages": {}}
    def _close(self, now, upto=None):
        for k, v in self.state["stages"].items():
            if k != upto and "end" not in v: v["end"] = now; v["status"] = "done"
    def set(self, stage, done=None, total=None, msg="", state=None):
        now = time.time(); s = self.state["stages"].get(stage)
        if s is None or "end" in s:                      # entering (or re-entering) a stage closes whatever was running
            self._close(now, upto=stage); s = self.state["stages"].setdefault(stage, {"start": now}); s.pop("end", None)
        s["status"] = "running"
        if done is not None: s["done"] = done
        if total is not None: s["total"] = total
        if msg: s["message"] = msg
        self.state.update(current=stage, updated=time.strftime("%Y-%m-%d %H:%M:%S"))
        if state: self.finish(state)
        else: self.p.write_text(json.dumps(self.state, indent=1))
        if msg: self.note(f"{stage}: {msg}")
    def note(self, msg):
        with open(self.log, "a") as f: f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    def finish(self, state):
        now = time.time(); self.state["state"] = state; self.state["ended"] = now
        for v in self.state["stages"].values():
            if "end" not in v: v["end"] = now; v["status"] = "done" if state == "finished" else "failed"
        self.p.write_text(json.dumps(self.state, indent=1))


def _sha(path) -> str:
    h = hashlib.sha256(); h.update(Path(path).read_bytes()); return h.hexdigest()[:16]


def prepare_receptor(src: Path, dst: Path) -> dict:
    """PDB/mmCIF -> protein-only PDB for docking. Waters, ligands, ions and cofactors are removed (documented limitation)."""
    from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
    parser = MMCIFParser(QUIET=True) if src.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    try:
        st = parser.get_structure("rec", str(src))
        n_chain_before = sum(1 for _ in st[0])
    except Exception as e:
        raise ValueError(f"Could not read {src.name} as a PDB/mmCIF protein structure ({type(e).__name__}: {e})") from e
    if n_chain_before == 0: raise ValueError(f"{src.name} contains no chains")
    class Prot(Select):
        def accept_residue(self, r): return r.id[0] == " "
        def accept_atom(self, a): return a.get_altloc() in (" ", "A")
    io = PDBIO(); io.set_structure(st); io.save(str(dst), Prot())
    n_res = sum(1 for l in open(dst) if l.startswith("ATOM") and l[12:16].strip() == "CA")
    het = sum(1 for r in st.get_residues() if r.id[0] not in (" ",) and r.id[0] != "W")
    if n_res == 0: raise ValueError("No protein residues found in the receptor file")
    longchain = [c.id for c in st[0] if len(c.id) > 1]
    if longchain: raise ValueError(f"Multi-character chain IDs {longchain} cannot be written to PDB for GNINA; renumber chains first")
    return {"n_chains": n_chain_before, "n_residues": n_res, "n_hetero_residues_removed": het}


def _ligand_to_pdb(path: Path, dst: Path) -> Path:
    if path.suffix.lower() == ".pdb": shutil.copy(path, dst); return dst
    m = (Chem.MolFromMolFile(str(path), removeHs=False) if path.suffix.lower() in (".sdf", ".mol")
         else Chem.MolFromMol2File(str(path), removeHs=False))
    if m is None: raise ValueError(f"Cannot read known ligand {path}")
    Chem.MolToPDBFile(m, str(dst)); return dst


def make_pockets(cfg: RunConfig, receptor: Path, rundir: Path, prog: Progress) -> tuple[list[Pocket], dict]:
    from . import pocket_tools
    pockets, status = [], {}
    if cfg.known_ligand:
        pdb = _ligand_to_pdb(Path(cfg.known_ligand), rundir / "inputs" / "known_ligand.pdb")
        pockets.append(from_known_ligand("known_ligand", str(receptor), str(pdb), padding=cfg.box_padding)); status["known_ligand"] = "ok"
    if cfg.box_center:
        pockets.append(from_user_box("user_box", str(receptor), cfg.box_center, tuple(cfg.box_size))); status["user_box"] = "ok"
    if cfg.predict_pockets or not pockets:
        for tool in cfg.predictors:
            prog.set(f"pocket:{tool}", msg=f"running {tool}")
            try:
                od = rundir / "pockets" / tool
                if tool == "p2rank": pk = pocket_tools.run_p2rank(str(receptor), od, cfg.max_pockets_per_tool, padding=4.0)
                elif tool == "pocketminer": pk, _ = pocket_tools.run_pocketminer(str(receptor), od, threshold=cfg.pocketminer_threshold, max_pockets=cfg.max_pockets_per_tool)
                elif tool == "lacuna": pk = pocket_tools.run_lacuna(str(receptor), od, top=cfg.max_pockets_per_tool)
                else: raise ValueError(f"unknown pocket tool {tool}")
                pockets += pk; status[tool] = f"ok ({len(pk)} pockets)" if pk else "ran, no pockets above threshold"
                prog.set(f"pocket:{tool}", done=1, total=1, msg=status[tool])
            except Exception as e:                       # documented: no silent substitution
                status[tool] = f"FAILED: {e}"; prog.set(f"pocket:{tool}", msg=status[tool])
    if not pockets:
        raise RuntimeError("No pockets available. Pocket tool status: " + json.dumps(status) + ". Supply a known ligand or box coordinates.")
    return pockets, status


def _covalent_stage(cfg: RunConfig, prepared: pd.DataFrame, pockets: list[Pocket], receptor: Path, rundir: Path, prog: Progress):
    """Create adduct-form ligands for compounds whose recognised warhead can react with a target residue type; dock them covalently.
    Every requested residue (cfg.residues()) is docked as its own covalent job set; rows are keyed by residue so they are never mixed."""
    ok = prepared[prepared.status == "ok"]
    cov_dir = rundir / "ligands_cov"; cov_dir.mkdir(exist_ok=True)
    forms_rows, notes = [], []
    for _, r in ok.iterrows():
        mol = Chem.RemoveHs(Chem.SDMolSupplier(r.ligand_sdf, removeHs=False)[0])
        forms, unsup = covalent_candidates(mol)
        for i, f in enumerate(forms):
            wh = get_warhead(f.warhead)
            m3, note = embed_3d(f.product)
            if m3 is None: notes.append({"ligand_key": r.ligand_key, "note": f"adduct 3D embedding failed ({f.warhead})"}); continue
            m3.SetProp("_Name", r.compound_id)
            p = cov_dir / f"{r.ligand_key}__cov{i}.sdf"
            with Chem.SDWriter(str(p)) as w: w.write(m3)
            forms_rows.append({"base": f"{r.ligand_key}__cov{i}", "ligand_key": r.ligand_key, "compound_id": r.compound_id, "warhead": f.warhead, "attach_idx": f.attach_idx,
                               "adduct_smiles": Chem.MolToSmiles(f.product), "adduct_sdf": str(p), "gnina_smarts": f.gnina_smarts, "typical_residues": ",".join(wh.typical_residues)})
        for u in unsup: notes.append({"ligand_key": r.ligand_key, "note": f"recognised warhead '{u}' has no covalent transform; noncovalent docking only"})
    rows, specs = [], {}
    for k, residue in enumerate(cfg.residues()):
        chain, num, _ = parse_residue(residue)
        for pk in pockets:
            rp = pk.receptor_path or str(receptor)
            try: atoms = residue_atoms(load_structure(rp), residue)
            except ValueError as e: notes.append({"ligand_key": "*", "note": f"pocket {pk.pocket_id}: {e}; covalent docking skipped for this pocket"}); continue
            resname = next(iter(atoms.values()))[1]
            names = NUCLEOPHILE_ATOMS.get(resname)
            if not names: notes.append({"ligand_key": "*", "note": f"target residue type {resname} has no defined nucleophile atom"}); continue
            nuc = next((n for n in names if n in atoms), None)
            if nuc is None: continue
            # covalent docking is only attempted where the nucleophile is near the pocket centre (search boxes can be 32 A wide and covalent jobs are slow);
            # elsewhere the job is skipped and recorded in covalent_report.csv
            dctr = float(np.linalg.norm(np.asarray(atoms[nuc][0], dtype=float) - np.asarray(pk.center)))
            if cfg.covalent_max_center_A and dctr > cfg.covalent_max_center_A:
                notes.append({"ligand_key": "*", "note": f"{residue} {nuc} is {dctr:.1f} A from the centre of {pk.pocket_id} (> {cfg.covalent_max_center_A:.0f} A, user-set speed limit); covalent docking not attempted there"}); continue
            for f in forms_rows:
                key = f["base"] + (f"__r{k}" if k else "")
                if f["typical_residues"].lower().find(resname.lower()[:3]) < 0:
                    notes.append({"ligand_key": f["ligand_key"], "note": f"{f['warhead']} is not typically reactive toward {resname}; covalent docking not attempted at {pk.pocket_id}"}); continue
                specs[(key, pk.pocket_id)] = {"rec_atom": f"{chain}:{num}:{nuc}", "lig_smarts": f["gnina_smarts"], "bond_order": 1}
                if not any(x["cov_key"] == key for x in rows):
                    rows.append({"cov_key": key, **{c: f[c] for c in ("ligand_key", "compound_id", "warhead", "attach_idx", "adduct_smiles", "adduct_sdf", "gnina_smarts", "typical_residues")}, "residue": residue})
    cmap = pd.DataFrame(rows)
    pd.DataFrame(notes).to_csv(rundir / "covalent_report.csv", index=False)
    if cmap.empty or not specs: return cmap
    cmap.to_csv(rundir / "covalent_map.csv", index=False)
    covdf = cmap.assign(status="ok", ligand_sdf=cmap.adduct_sdf)[["compound_id", "cov_key", "status", "ligand_sdf"]].rename(columns={"cov_key": "ligand_key"})
    prog.set("covalent_docking", 0, len({k[0] for k in specs}), f"covalent docking at {', '.join(cfg.residues())}")
    run_screen(covdf, str(receptor), pockets, cfg.dock_cfg(), rundir, mode=MODE_COVALENT, covalent_specs=specs, target_residue=cfg.target_residue,
               progress=lambda n, t, c, s: prog.set("covalent_docking", n, t, f"{c}: {s}"), pairs=set(specs))
    return cmap


def run(rundir: str | Path) -> None:
    rundir = Path(rundir); cfg = RunConfig.load(rundir)
    slot = None
    if hosting.HOSTED:                                  # public server: wait for a free machine-wide slot (state shows as Queued meanwhile)
        slot = hosting.Slot(); slot.acquire(); (rundir / "queued.json").unlink(missing_ok=True)
    prog = Progress(rundir, plan(cfg))
    try:
        (rundir / "inputs").mkdir(exist_ok=True)
        # 1. inputs
        prog.set("inputs", msg="validating receptor and library")
        src = Path(cfg.receptor); rec = rundir / "inputs" / "receptor.pdb"
        rinfo = prepare_receptor(src, rec)
        if cfg.residues():
            from .residues import validate_target
            specs_ok = []
            for r in cfg.residues():
                v = validate_target(str(rec), r)
                if not v["ok"]: raise ValueError(f"Target residue '{r}' is not usable: " + " ".join(v["errors"]))
                specs_ok.append(v["spec"])
            specs_ok = list(dict.fromkeys(specs_ok)); cfg.target_residue, cfg.target_residues = specs_ok[0], specs_ok[1:]
        lib = rundir / "inputs" / Path(cfg.library).name
        if not lib.exists(): shutil.copy(cfg.library, lib)
        manifest = {"pscreen": __version__, "python": sys.version.split()[0], "platform": platform.platform(), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
                    "receptor_sha256_16": _sha(src), "library_sha256_16": _sha(lib), "receptor_info": rinfo, "gnina": gnina_version(), "config": asdict(cfg)}
        # 2. compounds
        prog.set("prepare", msg="standardizing compounds")
        prepared = prepare_library(lib, rundir)
        summ = prep_summary(prepared); manifest["library_summary"] = summ
        prog.set("prepare", done=summ["total"], total=summ["total"], msg=f"library: {summ}")
        # 3. pockets
        prog.set("pockets", msg="defining pockets")
        pockets, pstat = make_pockets(cfg, rec, rundir, prog); manifest["pocket_tools"] = pstat
        (rundir / "pockets.json").write_text(json.dumps([asdict(p) for p in pockets], indent=1))
        prog.note(f"pockets: {len(pockets)} defined: {pstat}")
        # 4. noncovalent docking
        tr = cfg.target_residue or None
        n_jobs = int(((prepared.status == "ok").sum()) * len(pockets))
        prog.set("docking", 0, n_jobs, "starting GNINA docking")
        store = run_screen(prepared, str(rec), pockets, cfg.dock_cfg(), rundir, target_residue=tr,
                           progress=lambda n, t, c, s: prog.set("docking", n, t, f"{c}: {s}"))
        # 5. covalent
        cmap = pd.DataFrame()
        if cfg.covalent and cfg.residues():
            prog.set("covalent_docking", msg="building adduct-form ligands")
            try: cmap = _covalent_stage(cfg, prepared, pockets, rec, rundir, prog)
            except Exception as e:
                manifest["covalent_error"] = f"{e}"; prog.set("covalent_docking", msg=f"covalent stage failed (noncovalent results kept): {e}")
        # 6. ADMET (optional)
        props = None
        if cfg.chemprop:
            prog.set("properties", msg="ADMET-AI predictions")
            props, meta = chem_props.predict(prepared, rundir / "props", cfg.chemprop_endpoints); manifest["chemprop"] = meta
            if props is None: prog.set("properties", msg=meta.get("disabled", "disabled"))
        # 7. results + exports
        prog.set("results", msg="assembling results")
        if store.db.execute("SELECT COUNT(*) FROM jobs WHERE status='done'").fetchone()[0] == 0:
            _failed = [r[0] for r in store.db.execute("SELECT error FROM jobs WHERE status='failed' LIMIT 1")]
            raise RuntimeError("No compound could be docked, so there are no results. First docking error: " + (_failed[0] if _failed else "none recorded") + ". Check the protein, the compound structures and the binding site.")
        exp = pd.read_csv(cfg.experimental_csv) if cfg.experimental_csv else None
        res = build_results(prepared, store, pockets, str(rec), tr, exp, cmap if len(cmap) else None, props)
        res.to_csv(rundir / "results_all_poses.csv", index=False)
        best_pose_per_compound(res, "gnina_cnn_vs", None).to_csv(rundir / "results_best_pose.csv", index=False)
        export_structures(res, rundir)
        jobs = store.jobs(); manifest["jobs"] = jobs.groupby(["mode", "status"]).size().rename("n").reset_index().to_dict("records")
        manifest["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
        (rundir / "run_manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
        try:
            from .api import export_outputs
            export_outputs(rundir)
        except Exception as e: prog.note(f"output package not written: {e}")
        prog.set("results", msg="finished", state="finished")
    except Exception as e:
        (rundir / "error.txt").write_text(traceback.format_exc())
        prog.set(prog.state.get("current", "run"), msg=f"FAILED: {e}", state="failed")
        raise
    finally:
        if slot: slot.release()


def export_structures(res: pd.DataFrame, rundir: Path) -> Path:
    """Bundle: all pose SDFs referenced in the results table, receptor conformations, pockets.json, results CSVs."""
    z = rundir / "structures_bundle.zip"
    docked = res[res.processing_status == "docked"] if "processing_status" in res else res
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(set(docked.pose_sdf.dropna())): zf.write(p, f"poses/{Path(p).name}")
        for f in ["inputs/receptor.pdb", "pockets.json", "results_best_pose.csv", "run_config.json", "run_manifest.json"]:
            if (rundir / f).exists(): zf.write(rundir / f, f)
        for f in (rundir / "pockets").glob("lacuna/conformer_*.pdb"): zf.write(f, f"receptor_conformers/{f.name}")
    return z


if __name__ == "__main__":
    run(sys.argv[1])
