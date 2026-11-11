"""Hosted (public web) mode: private per-visitor workspaces, free-tier caps, a global run queue and automatic clean-up.

Off by default, so the local app is unchanged. Enable with  PSCREEN_HOSTED=1  (the Dockerfile does).

- Workspace: a random 128-bit token kept in the page URL (?ws=...). It is the only thing that grants access to a visitor's runs, so a
  bookmarked URL returns to the same experiments and nobody else can list them. No accounts, no cookies, no tracking.
- Caps: server-side limits on library size, protein size, residues, pockets, settings and concurrent runs per workspace (`check_config`).
- Queue: every pipeline process first takes one of N machine-wide slots (flock files). Waiting runs show as "queued" with their place.
- Clean-up: workspaces untouched for PSCREEN_RETENTION_H hours are deleted (`purge_old`).
"""
from __future__ import annotations
import fcntl, json, os, re, secrets, shutil, sqlite3, time
from pathlib import Path

HOSTED = os.environ.get("PSCREEN_HOSTED") == "1"
DATA = Path(os.environ.get("PSCREEN_DATA", "/tmp/pscreen_data" if HOSTED else "."))
UNLIMITED = HOSTED and os.environ.get("PSCREEN_LIMITS", "on") == "off"   # operator's own machine: no caps, runs behave as in the local app
RETENTION_H = float(os.environ.get("PSCREEN_RETENTION_H", "48"))
SLOTS = int(os.environ.get("PSCREEN_SLOTS", "1"))            # pipelines that may run at once on this machine
EXAMPLES = Path(__file__).resolve().parent.parent / "examples"   # read-only demo runs copied into every new workspace

LIMITS = {
    "max_compounds": int(os.environ.get("PSCREEN_MAX_COMPOUNDS", "25")),
    "max_protein_mb": 15,
    "max_library_mb": 5,
    "max_residues": 2,                # covalent target residues
    "max_pockets_per_tool": 2,
    "max_exhaustiveness": 8,
    "max_active_runs_per_workspace": 1,
    "max_runs_per_workspace": 6,
    "timeout_s_per_job": 600,
    "cpu_per_job": 2,
    "workers": 1,
}
_TOKEN = re.compile(r"^[0-9a-f]{32}$")


def safe_name(name: str) -> str:
    """A client-supplied upload name reduced to a harmless file name (no folders, no '..', limited characters, keeps the extension)."""
    base = Path(str(name).replace("\\", "/")).name
    base = re.sub(r"[^A-Za-z0-9._\-]", "_", base).lstrip(".")[:80]
    return base or "upload"


MAX_QUEUE = int(os.environ.get("PSCREEN_MAX_QUEUE", "6"))      # pipelines waiting or running across ALL workspaces


def global_queue_len() -> int:
    if not HOSTED: return 0
    return sum(1 for q in (DATA / "workspaces").glob("*/runs/*/pid") if _pid_alive(q))


def new_token() -> str:
    return secrets.token_hex(16)


def valid_token(t: str | None) -> bool:
    return bool(t) and bool(_TOKEN.match(t))


def workspace_dir(token: str) -> Path:
    return DATA / "workspaces" / token


def runs_root(token: str | None) -> Path:
    """Where this visitor's runs live. Local mode keeps the original `runs/` folder."""
    if not HOSTED: return Path("runs")
    if not valid_token(token): raise ValueError("invalid workspace token")
    return workspace_dir(token) / "runs"


def ensure_workspace(token: str) -> Path:
    """Create the workspace (first visit) and seed it with the example run(s). Touches the folder so clean-up sees it as active."""
    root = runs_root(token); first = not root.exists(); root.mkdir(parents=True, exist_ok=True)
    if first and EXAMPLES.exists():
        for ex in EXAMPLES.iterdir():
            if ex.is_dir() and (ex / "run_config.json").exists(): shutil.copytree(ex, root / ex.name, dirs_exist_ok=True); relocate(root / ex.name)
    touch(token)
    return root


PLACEHOLDER = "@RUN@"          # stands for the run folder inside shipped example runs (no machine-specific paths in the repo)


def relocate(rundir: Path):
    """Replace the @RUN@ placeholder in a copied example run by its real folder (text files and the results database)."""
    here = str(rundir.resolve())
    for f in rundir.rglob("*"):
        if f.is_file() and f.suffix in (".json", ".csv", ".txt", ".md", ".log") and f.stat().st_size < 50e6:
            t = f.read_text(errors="ignore")
            if PLACEHOLDER in t: f.write_text(t.replace(PLACEHOLDER, here))
    db = rundir / "results.sqlite"
    if db.exists():
        c = sqlite3.connect(db)
        c.execute("UPDATE poses SET pose_sdf = replace(pose_sdf, ?, ?)", (PLACEHOLDER, here))
        c.execute("UPDATE meta SET value = replace(value, ?, ?)", (PLACEHOLDER, here)); c.commit(); c.close()


def touch(token: str):
    if HOSTED: (workspace_dir(token) / ".active").write_text(str(time.time()))


def purge_old(now: float | None = None) -> list[str]:
    """Delete workspaces idle longer than RETENTION_H that have no running pipeline. Returns the removed tokens."""
    now = now or time.time(); gone = []
    base = DATA / "workspaces"
    if not HOSTED or not base.exists(): return gone
    for w in base.iterdir():
        if not (w.is_dir() and valid_token(w.name)): continue
        try: last = float((w / ".active").read_text())
        except (OSError, ValueError): last = w.stat().st_mtime
        busy = any(_pid_alive(p) for p in w.glob("runs/*/pid"))
        if not busy and now - last > RETENTION_H * 3600: shutil.rmtree(w, ignore_errors=True); gone.append(w.name)
    return gone


def _pid_alive(pidfile: Path) -> bool:
    try: os.kill(int(pidfile.read_text()), 0); return True
    except (OSError, ValueError): return False


# ------------------------------------------------------------------------------------------------------------------ caps
def check_config(cfg_kw: dict, n_compounds: int | None, protein_bytes: int, library_bytes: int, active_runs: int, total_runs: int) -> list[str]:
    """Returns human-readable reasons a run is refused (empty = allowed). Only enforced in hosted mode."""
    if not HOSTED or UNLIMITED: return []
    L, e = LIMITS, []
    if n_compounds and n_compounds > L["max_compounds"]: e.append(f"This free online version docks at most {L['max_compounds']} compounds per run (yours has {n_compounds}). Install pscreen locally for larger libraries.")
    if protein_bytes > L["max_protein_mb"] * 2**20: e.append(f"Protein file is larger than {L['max_protein_mb']} MB.")
    if library_bytes > L["max_library_mb"] * 2**20: e.append(f"Compound file is larger than {L['max_library_mb']} MB.")
    if len(cfg_kw.get("target_residues", [])) + bool(cfg_kw.get("target_residue")) > L["max_residues"]: e.append(f"At most {L['max_residues']} covalent target residues per run online.")
    if active_runs >= L["max_active_runs_per_workspace"]: e.append("You already have a run in progress or queued. Wait for it to finish (you can close the tab and return with this page's URL).")
    if global_queue_len() >= MAX_QUEUE: e.append("The server is busy (too many runs queued). Please try again later.")
    if total_runs >= L["max_runs_per_workspace"]: e.append(f"This workspace holds {L['max_runs_per_workspace']} experiments already. Delete one or open a new workspace.")
    return e


def clamp_config(cfg_kw: dict) -> dict:
    """Apply the free-tier ceilings to compute settings (never raises them)."""
    if not HOSTED or UNLIMITED: return cfg_kw
    L = LIMITS; k = dict(cfg_kw)
    k["max_pockets_per_tool"] = min(int(k.get("max_pockets_per_tool", 3)), L["max_pockets_per_tool"])
    k["exhaustiveness"] = min(int(k.get("exhaustiveness", 8)), L["max_exhaustiveness"])
    k["cpu_per_job"] = min(int(k.get("cpu_per_job", 3)), L["cpu_per_job"])
    k["workers"] = min(int(k.get("workers", 2)), L["workers"])
    k["timeout_s"] = min(int(k.get("timeout_s", 1800)), L["timeout_s_per_job"])
    return k


# ------------------------------------------------------------------------------------------------------------------ queue
class Slot:
    """One of SLOTS machine-wide run slots, held for the life of the pipeline process (the OS releases it if the process dies)."""
    def __init__(self): self.f = None

    def acquire(self, on_wait=None, poll_s: float = 2.0):
        if UNLIMITED: return                                  # no queue: every run starts immediately
        d = DATA / "slots"; d.mkdir(parents=True, exist_ok=True)
        while True:
            for i in range(SLOTS):
                f = open(d / f"slot{i}.lock", "a+")
                try: fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB); self.f = f; return
                except OSError: f.close()
            if on_wait: on_wait()
            time.sleep(poll_s)

    def release(self):
        if self.f: fcntl.flock(self.f, fcntl.LOCK_UN); self.f.close(); self.f = None


def queue_position(rundir: Path) -> int:
    """1 = next. Counts live pipelines of other workspaces that were queued earlier (by their queued.json time)."""
    me = _queued_at(rundir)
    if me is None: return 0
    ahead = 0
    for q in (DATA / "workspaces").glob("*/runs/*/queued.json"):
        if q.parent == rundir: continue
        t = _queued_at(q.parent)
        if t is not None and t < me and _pid_alive(q.parent / "pid"): ahead += 1
    return ahead + 1


def _queued_at(rundir: Path) -> float | None:
    try: return float(json.loads((rundir / "queued.json").read_text())["t"])
    except (OSError, ValueError, KeyError): return None
