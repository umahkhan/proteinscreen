"""Workflow builder: a node canvas (palette | canvas) beside an inputs panel, with a live per-step log and timers.

The canvas is a small bidirectional Streamlit component (wf_component/index.html). Its value is the list of nodes on the canvas; `to_config`
turns that graph into a RunConfig, so every combination of pocket models / docking lanes / property models the user builds is exactly what runs.
"""
from __future__ import annotations
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components
from .. import chem_props, pocket_tools
from ..pipeline import RunConfig, plan, prepare_receptor, stage_label
from .. import hosting
from ..libformat import FORMAT_HELP, TEMPLATE_CSV, validate_library
from ..runio import list_runs as _list_runs
from .layout import viewport_h



def runs_root() -> Path:
    """This visitor's runs folder (hosted: private workspace from the ?ws= token; local: ./runs)."""
    return hosting.runs_root(st.session_state.get("ws"))


def list_runs() -> list[Path]:
    return _list_runs(runs_root())

_comp = components.declare_component("pscreen_workflow", path=str(Path(__file__).parent / "wf_component"))

BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
NODE_DEFS = [
    dict(type="protein", group=1, title="Protein structure", desc="PDB / mmCIF (Inputs panel)", color="#374151", locked=True, initial=True, singleton=True, params=[], blurb="Uploaded in the Inputs panel →"),
    dict(type="library", group=1, title="Compound library", desc="SMILES / CSV / SDF (Inputs panel)", color="#374151", locked=True, initial=True, singleton=True, params=[], blurb="Uploaded in the Inputs panel →"),
    dict(type="p2rank", group=2, title="P2Rank", desc="Conventional pockets (geometry + ML)", color=ORANGE, singleton=True,
         params=[dict(k="max_pockets", label="Max pockets", v=3, min=1)], blurb="Ligandable surface pockets."),
    dict(type="pocketminer", group=2, title="PocketMiner", desc="Cryptic-pocket propensity per residue", color=GREEN, singleton=True,
         params=[dict(k="threshold", label="Prob. ≥", v=0.7, step=0.05, min=0.3), dict(k="max_pockets", label="Max pockets", v=3, min=1)], blurb="Residue GNN, clustered into pockets (our threshold)."),
    dict(type="lacuna", group=2, title="Lacuna", desc="Conformational ensemble pockets", color=BLUE, singleton=True,
         params=[dict(k="max_pockets", label="Max pockets", v=3, min=1)], blurb="NMA ensemble; docks into the most open conformer."),
    dict(type="known_ligand", group=2, title="Known-ligand site", desc="Site from a co-crystal ligand file", color="#6250d6", singleton=True, params=[], blurb="Uses the known-ligand file in the Inputs panel."),
    dict(type="box", group=2, title="Box coordinates", desc="Manual search box", color="#6250d6", singleton=True,
         params=[dict(k="cx", label="Centre x", v=0.0, step=0.5), dict(k="cy", label="Centre y", v=0.0, step=0.5), dict(k="cz", label="Centre z", v=0.0, step=0.5), dict(k="size", label="Size Å", v=22.0, step=1, min=8)], blurb=""),
    dict(type="target", group=3, title="Target residue(s)", desc="Residues for covalent docking + geometry", color="#b4232a", singleton=True,
         params=[dict(k="residues", label="Residues", type="text", v="", ph="A:CYS727, A:CYS106")], blurb="One or more, comma-separated (CHAIN:RESNAME+NUM)."),
    dict(type="gnina", group=4, title="GNINA docking", desc="Non-covalent docking (required)", color="#1d4ed8", singleton=True, initial=True,
         params=[dict(k="exhaustiveness", label="Exhaust.", v=8, min=1), dict(k="num_modes", label="Poses", v=9, min=1), dict(k="cpu", label="CPUs/job", v=3, min=1), dict(k="workers", label="Parallel", v=2, min=1)], blurb=""),
    dict(type="gnina_cov", group=4, title="GNINA covalent", desc="Adduct-form covalent docking", color="#c2410c", singleton=True, needs=["target"],
         params=[dict(k="max_center", label="Limit Å (0=off)", v=0, min=0)], blurb="Speed limit skips pockets farther than this from the residue (biases toward that site)."),
    dict(type="warhead", group=5, title="Warhead + RDKit properties", desc="Always on", color="#4b5563", locked=True, initial=True, singleton=True, params=[], blurb="Warhead classes, MW, cLogP, TPSA…"),
    dict(type="admet", group=5, title="ADMET-AI", desc="Chemprop ADMET predictions", color="#0f766e", singleton=True, params=[], blurb="Property predictions, not binding or reactivity."),
    dict(type="results", group=6, title="Results", desc="Ranked table + 3D explorer", color="#111827", locked=True, initial=True, singleton=True, params=[], blurb="Opens in the Results tab."),
]
PRESETS = {
    "Blank (docking only)": [dict(type="gnina")],
    "P2Rank + docking (no ADMET)": [dict(type="p2rank"), dict(type="gnina")],
    "P2Rank + PocketMiner + docking + ADMET": [dict(type="p2rank"), dict(type="pocketminer"), dict(type="gnina"), dict(type="admet")],
    "Covalent: P2Rank pocket + covalent docking (quick, no ADMET)": [dict(type="p2rank"), dict(type="target"), dict(type="gnina"), dict(type="gnina_cov")],
    "Covalent: cryptic pockets (PocketMiner + Lacuna) + ADMET": [dict(type="pocketminer"), dict(type="lacuna"), dict(type="target"), dict(type="gnina"), dict(type="gnina_cov"), dict(type="admet")],
    "Covalent: all three pocket models + ADMET": [dict(type="p2rank"), dict(type="pocketminer"), dict(type="lacuna"), dict(type="target"), dict(type="gnina"), dict(type="gnina_cov"), dict(type="admet")],
}
STAGE_NODES = {"inputs": ["protein", "library"], "prepare": ["warhead"], "pockets": ["known_ligand", "box"], "pocket:p2rank": ["p2rank"], "pocket:pocketminer": ["pocketminer"],
               "pocket:lacuna": ["lacuna"], "docking": ["gnina"], "covalent_docking": ["gnina_cov"], "properties": ["admet"], "results": ["results"]}


# ---------------------------------------------------------------------------------------------------------------- run control
def pid_alive(pid: int) -> bool:
    try: os.kill(pid, 0); return True
    except OSError: return False


def launch(rundir: Path):
    if hosting.HOSTED: (rundir / "queued.json").write_text(json.dumps({"t": time.time()}))
    log = open(rundir / "stdout.log", "ab")
    p = subprocess.Popen([sys.executable, "-m", "pscreen.pipeline", str(rundir)], stdout=log, stderr=log, start_new_session=True)
    (rundir / "pid").write_text(str(p.pid))


def stop(rundir: Path):
    pf = rundir / "pid"
    if pf.exists():
        try: os.killpg(int(pf.read_text()), signal.SIGTERM)
        except (OSError, ValueError): pass
    if sys.platform == "darwin":                        # GNINA runs in the Lima VM on macOS; match the run's full path so other visitors' runs of the same name are untouched
        subprocess.run(["limactl", "shell", "gnina", "--", "pkill", "-f", f"{rundir.resolve()}/"], capture_output=True, timeout=30)


ACTIVE = ("running", "queued")        # states that keep the page refreshing and count against the per-visitor active-run cap


def queue_note(rd: Path) -> str:
    n = hosting.queue_position(rd) if hosting.HOSTED else 0
    return f"Waiting for a free server slot (place {n} in the queue). You can close this tab and come back with the same page URL." if n else ""


def run_state(rd: Path) -> tuple[str, dict | None]:
    """('running'|'finished'|'failed'|'interrupted'|'new', progress dict)"""
    pf = rd / "progress.json"
    pid = int((rd / "pid").read_text()) if (rd / "pid").exists() else None
    alive = bool(pid and pid_alive(pid))
    if not pf.exists(): return (("queued" if (rd / "queued.json").exists() else "running") if alive else "new"), None
    try: d = json.loads(pf.read_text())
    except json.JSONDecodeError: return "running", None
    s = d.get("state", "running")
    if s == "running" and not alive: s = "interrupted"
    return s, d


# ---------------------------------------------------------------------------------------------------------------- graph -> config
def parse_residues(text: str) -> list[str]:
    return [t for t in re.split(r"[,;\s]+", text or "") if t]


def to_config(nodes: list[dict], name: str, paths: dict, receptor_pdb: str | None) -> tuple[dict | None, list[str], list[str]]:
    """Returns (RunConfig kwargs or None, errors, warnings)."""
    by = {n["type"]: n for n in nodes}; err, warn = [], []
    P = lambda t, k, d=None: by[t]["params"].get(k, d) if t in by else d
    if not paths.get("receptor"): err.append("Upload a protein structure (Inputs panel).")
    if not paths.get("library"): err.append("Upload a compound library (Inputs panel).")
    if "gnina" not in by: err.append("Add the GNINA docking step — without it there is nothing to rank.")
    sites = [t for t in ("p2rank", "pocketminer", "lacuna", "known_ligand", "box") if t in by]
    if not sites: err.append("Add at least one binding-site step (P2Rank, PocketMiner, Lacuna, known ligand or box).")
    if "known_ligand" in by and not paths.get("known_ligand"): err.append("The Known-ligand site step needs a ligand file (Inputs panel).")
    residues = parse_residues(P("target", "residues", "")) if "target" in by else []
    if "target" in by and not residues: warn.append("Target residue step is empty: no covalent docking or reactive-geometry columns.")
    if "gnina_cov" in by and not residues: err.append("Covalent docking needs at least one target residue.")
    if "admet" in by and not chem_props.available(): err.append("ADMET-AI is not installed (docs/SETUP.md).")
    for t, ok in (("pocketminer", pocket_tools.pm_possible()), ("lacuna", pocket_tools.lac_python() is not None)):
        if t in by and not ok: err.append(f"{t} is not installed (docs/SETUP.md).")
    if not name: err.append("Give the run a name.")
    elif (runs_root() / name / "run_config.json").exists(): err.append(f"A run named '{name}' already exists.")
    if err: return None, err, warn
    pred = [t for t in ("p2rank", "pocketminer", "lacuna") if t in by]
    kw = dict(name=name, receptor=paths["receptor"], library=paths["library"], known_ligand=paths.get("known_ligand", "") if "known_ligand" in by else "",
              experimental_csv=paths.get("experimental_csv", ""), target_residue=residues[0] if residues else "", target_residues=residues[1:],
              box_center=[P("box", "cx", 0.0), P("box", "cy", 0.0), P("box", "cz", 0.0)] if "box" in by else [], box_size=[float(P("box", "size", 22.0))] * 3,
              predict_pockets=bool(pred), predictors=pred, max_pockets_per_tool=int(max([P(t, "max_pockets", 3) for t in pred] or [3])),
              pocketminer_threshold=float(P("pocketminer", "threshold", 0.7)), exhaustiveness=int(P("gnina", "exhaustiveness", 8)), num_modes=int(P("gnina", "num_modes", 9)),
              cpu_per_job=int(P("gnina", "cpu", 3)), workers=int(P("gnina", "workers", 2)), covalent=("gnina_cov" in by and bool(residues)),
              covalent_max_center_A=float(P("gnina_cov", "max_center", 0) or 0), chemprop="admet" in by)
    return kw, err, warn


# ---------------------------------------------------------------------------------------------------------------- UI pieces
def node_status(rd: Path | None) -> dict:
    if rd is None: return {}
    state, d = run_state(rd)
    if not d: return {}
    out = {}
    for k, v in d["stages"].items():
        for t in STAGE_NODES.get(k, []): out[t] = "failed" if (state in ("failed", "interrupted") and v.get("status") == "running") else v.get("status", "")
    for k in d.get("plan", []):                       # planned-but-not-started steps stay unmarked
        pass
    return out


def fmt_t(sec: float | None) -> str:
    if sec is None: return "—"
    m, s = divmod(sec, 60); h, m = divmod(int(m), 60)
    return f"{h}:{m:02d}:{s:04.1f}" if h else f"{m}:{s:04.1f}"


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")[:60]


# ---------------------------------------------------------------------------------------------------------------- run views
STATE_CHIP = {"running": ("#e8f0fb", "#1b5aa8", "Running"), "finished": ("#e6f5ee", "#0f6c4a", "Finished"), "failed": ("#fdeaea", "#a12c2b", "Failed"),
              "interrupted": ("#fff3df", "#8a5a00", "Interrupted"), "new": ("#eceff3", "#55544f", "Not started"), "queued": ("#eceff3", "#55544f", "Queued")}


def chip(state: str) -> str:
    bg, fg, label = STATE_CHIP.get(state, STATE_CHIP["new"])
    return f'<span class="ps-chip" style="background:{bg};color:{fg}">{label}</span>'


def run_title(rd: Path) -> str:
    try: return json.loads((rd / "run_config.json").read_text()).get("name") or rd.name
    except Exception: return rd.name


def run_summary(rd: Path) -> str:
    """One line describing what the run contains, so experiments can be told apart."""
    try: c = json.loads((rd / "run_config.json").read_text())
    except Exception: return ""
    bits = []
    if c.get("known_ligand"): bits.append("known site")
    if c.get("box_center"): bits.append("box")
    names = {"p2rank": "P2Rank", "pocketminer": "PocketMiner", "lacuna": "Lacuna"}
    if c.get("predict_pockets") or not bits: bits += [names.get(t, t) for t in c.get("predictors", [])]
    res = [r for r in dict.fromkeys([c.get("target_residue"), *c.get("target_residues", [])]) if r]
    short = lambda r: (r.split(":")[-1][:1] + r.split(":")[-1][1:].lower()) if ":" in r else r
    out = " + ".join(bits)
    if c.get("covalent") and res: out += " · covalent " + ", ".join(short(r) for r in res)
    elif res: out += " · target " + ", ".join(short(r) for r in res)
    if c.get("chemprop"): out += " · ADMET"
    return out


def run_nodes_html(rd: Path) -> str:
    """Condensed chain of the optional workflow nodes this experiment used, in pipeline order (colours match the canvas)."""
    nodes = []
    try: nodes = json.loads((rd / "workflow.json").read_text())["nodes"]
    except Exception:
        try:
            c = json.loads((rd / "run_config.json").read_text())
            nodes = [{"type": t, "params": {}} for t in c.get("predictors", [])] + ([{"type": "known_ligand", "params": {}}] if c.get("known_ligand") else []) + ([{"type": "box", "params": {}}] if c.get("box_center") else [])
            res = [r for r in [c.get("target_residue"), *c.get("target_residues", [])] if r]
            if res: nodes.append({"type": "target", "params": {"residues": ", ".join(res)}})
            nodes += [{"type": "gnina", "params": {}}] + ([{"type": "gnina_cov", "params": {}}] if c.get("covalent") else []) + ([{"type": "admet", "params": {}}] if c.get("chemprop") else [])
        except Exception: return ""
    defs = {d["type"]: d for d in NODE_DEFS}; out = []
    for n in sorted(nodes, key=lambda n: (defs.get(n["type"], {}).get("group", 9))):
        d = defs.get(n["type"])
        if not d or d.get("locked"): continue
        extra = f' · {n["params"].get("residues")}' if n["type"] == "target" and n["params"].get("residues") else ""
        out.append(f'<span class="ps-node" style="border-left-color:{d["color"]}">{d["title"]}{extra}</span>')
    return '<span class="ps-arrow">›</span>'.join(out)


def started_str(rd: Path, d: dict | None) -> str:
    t = (d or {}).get("started") or ((rd / "run_config.json").stat().st_mtime if (rd / "run_config.json").exists() else None)
    return time.strftime("%b %d · %H:%M", time.localtime(t)) if t else "—"


def go(page: str, run: str | None = None):
    """Navigate from inside a widget handler (the nav widget is already instantiated, so app.py applies this at the top of the next run)."""
    st.session_state["_goto"] = {"nav": page, **({"res_run": run} if run else {})}
    st.rerun()


def run_log_card(rd: Path):
    state, d = run_state(rd); now = time.time()
    if d is None:
        try: d = {"plan": plan(RunConfig.load(rd)), "stages": {}, "started": now}
        except Exception: d = {"plan": [], "stages": {}, "started": now}
    d.setdefault("started", min([v["start"] for v in d["stages"].values() if "start" in v] or [(rd / "progress.json").stat().st_mtime if (rd / "progress.json").exists() else now]))
    stages = d["stages"]
    end = d.get("ended") or (now if state == "running" else max([v.get("end", 0) for v in stages.values()] + [d["started"]]))
    total = end - d["started"] if (state == "running" or d.get("ended") or any("start" in v for v in stages.values())) else None
    st.markdown(f'''<div class="ps-card-h"><div><div class="ps-kick">Run log</div><div class="ps-big">{run_title(rd)}</div>
        <div class="ps-meta">{run_summary(rd)}</div><div class="ps-meta">Run ID <code>{rd.name}</code> · started {started_str(rd, d)}</div></div>
        <div class="ps-right">{chip(state)}<div class="ps-kick" style="margin-top:8px">Total time</div><div class="ps-timer">{fmt_t(total)}</div></div></div>''', unsafe_allow_html=True)
    if state == "queued" and queue_note(rd): st.info(queue_note(rd))
    rows = ""
    order = list(d.get("plan", [])) + [k for k in stages if k not in d.get("plan", [])]
    for i, k in enumerate(order, 1):
        v = stages.get(k); s = (v.get("status") or ("running" if state == "running" and k == d.get("current") else "done")) if v else "pending"
        if state in ("failed", "interrupted") and s == "running": s = "failed"
        el = ((v.get("end") or now) - v["start"]) if v and "start" in v else None
        bg, fg, lab = {"pending": ("#eceff3", "#55544f", "Pending"), "running": ("#e8f0fb", "#1b5aa8", "Running"), "done": ("#e6f5ee", "#0f6c4a", "Done"), "failed": ("#fdeaea", "#a12c2b", "Failed")}[s]
        done, tot = (v or {}).get("done"), (v or {}).get("total")
        if s == "done" and tot: done = tot            # a finished stage has completed all of its jobs (a resumed run skips jobs already done)
        pct = min(100, 100 * (done or 0) / max(1, tot)) if tot else (100 if s == "done" else 0)
        bar = (f'<div class="ps-bar"><div style="width:{pct}%;background:{"#2a78d6" if s == "running" else "#008300" if s == "done" else "#c9ced8"}"></div></div>'
               + (f'<span class="ps-cnt">{done or 0}/{tot}</span>' if tot else ""))
        msg = ((v or {}).get("message") or "")[:120].replace("<", "&lt;")
        rows += (f'<tr><td class="n">{i}</td><td class="st"><b>{stage_label(k)}</b></td><td><span class="ps-chip s" style="background:{bg};color:{fg}">{lab}</span></td><td>{bar}</td>'
                 f'<td class="tm">{fmt_t(el)}</td><td class="msg">{msg}</td></tr>')
    st.markdown(f'<table class="ps-steps"><thead><tr><th></th><th>Step</th><th>Status</th><th>Progress</th><th style="text-align:right">Time</th><th>Latest message</th></tr></thead><tbody>{rows}</tbody></table>', unsafe_allow_html=True)
    if (rd / "error.txt").exists(): st.error("Run failed"); st.code((rd / "error.txt").read_text()[-1200:], language=None)
    with st.expander("Raw log", expanded=False):
        if (rd / "run.log").exists(): st.code("".join((rd / "run.log").read_text().splitlines(True)[-16:]), language=None)
    b1, b2, _ = st.columns([1, 1.3, 5])
    if state == "running" and b1.button("■ Stop run", key=f"stop_{rd.name}"): stop(rd); st.rerun()
    if state in ("failed", "interrupted") and b1.button("↻ Resume run", key=f"res_{rd.name}"): (rd / "error.txt").unlink(missing_ok=True); launch(rd); st.rerun()
    if state == "finished" and (rd / "results.sqlite").exists() and b2.button("Open results →", type="primary", key=f"open_{rd.name}"): go("Results", rd.name)


def experiments_list(active: str | None):
    runs = list_runs()
    st.markdown('<div class="ps-kick" style="margin:.2rem 0 .4rem">Past experiments</div>', unsafe_allow_html=True)
    hdr = st.columns([2.6, 1.0, 2.0, 1.0, 1.3, 1.9], vertical_alignment="center")
    for c, t in zip(hdr, ["Experiment", "Status", "Current / last step", "Elapsed", "Started", ""]):
        c.markdown(f'<div class="ps-th">{t}</div>', unsafe_allow_html=True)
    now = time.time()
    for rd in runs[:30]:
        state, d = run_state(rd)
        cur = (d or {}).get("current"); v = ((d or {}).get("stages") or {}).get(cur or "", {})
        started, end = (d or {}).get("started"), (d or {}).get("ended")
        el = ((end or now) - started) if started and (end or state == "running") else None
        done, tot = v.get("done"), v.get("total")
        step = (stage_label(cur) + (f" · {done or 0}/{tot}" if tot and state == "running" else "")) if cur else "—"
        viewing = rd.name == active
        c = st.columns([2.6, 1.0, 2.0, 1.0, 1.3, 1.9], vertical_alignment="center")
        c[0].markdown(f'<div class="ps-row{" sel" if viewing else ""}"><b>{run_title(rd)}</b><div class="ps-meta">{run_summary(rd)}</div></div>', unsafe_allow_html=True)
        c[1].markdown(chip(state), unsafe_allow_html=True)
        c[2].markdown(f'<div class="ps-cell">{step}</div>', unsafe_allow_html=True)
        c[3].markdown(f'<div class="ps-cell tab">{fmt_t(el)}</div>', unsafe_allow_html=True)
        c[4].markdown(f'<div class="ps-cell">{started_str(rd, d)}</div>', unsafe_allow_html=True)
        a1, a2 = c[5].columns(2)
        if a1.button("Log", key=f"lg_{rd.name}", disabled=viewing, use_container_width=True): st.session_state["active_run"] = rd.name; st.rerun()
        if a2.button("Results", key=f"rs_{rd.name}", disabled=not (state == "finished" and (rd / "results.sqlite").exists()), use_container_width=True): go("Results", rd.name)


def runs_page():
    """Run log first (the experiment being watched), then every past experiment with status, step and elapsed time."""
    S = st.session_state
    runs = list_runs()
    if not runs:
        st.markdown('<div class="ps-empty"><b>No experiments yet</b><br>Build a workflow, name the experiment and press <i>Run analysis</i>. Its log and timers appear here.</div>', unsafe_allow_html=True); return
    names = [r.name for r in runs]
    if S.get("active_run") not in names: S["active_run"] = names[0]
    live = any(run_state(r)[0] in ACTIVE for r in runs)

    def body():
        active = S["active_run"]
        run_log_card(runs_root() / active)
        st.markdown('<div class="ps-div"></div>', unsafe_allow_html=True)
        experiments_list(active)
        if not any(run_state(r)[0] in ACTIVE for r in list_runs()) and S.get("_runs_live"):
            S["_runs_live"] = False; st.rerun()
    S["_runs_live"] = live
    (st.fragment(run_every=1)(body) if live else body)()


# ---------------------------------------------------------------------------------------------------------------- build page
def render():
    """Workflow page: [workflow nodes | canvas] [inputs] — pressing Run analysis hands over to the Runs page."""
    S = st.session_state
    S.setdefault("wf_token", 0); S.setdefault("wf_preset", None); S.setdefault("wf_nodes", []); S.setdefault("wf_pvals", {})
    RUNS = runs_root(); RUNS.mkdir(parents=True, exist_ok=True)
    active = RUNS / S["active_run"] if S.get("active_run") and (RUNS / S["active_run"]).exists() else None
    live_runs = [r for r in list_runs() if run_state(r)[0] in ACTIVE]
    vh = viewport_h()

    top = st.columns([4.0, 1.3], gap="medium")
    with top[0]:
        t1, t2 = st.columns([1.1, 2.4], vertical_alignment="center")
        t1.markdown('<div class="ps-kick">Workflow canvas</div>', unsafe_allow_html=True)
        def _preset():
            ch = S.get("wf_preset_pick")
            if ch in PRESETS: S["wf_token"] += 1; S["wf_preset"] = dict(token=S["wf_token"], nodes=PRESETS[ch])
        t2.selectbox("Start from", [None, *PRESETS], key="wf_preset_pick", on_change=_preset, label_visibility="collapsed",
                     format_func=lambda x: "Start from a preset…" if x is None else x, help="Presets only place nodes on the canvas — you can add, remove and edit nodes afterwards.")
        status = {}                      # node colouring follows the selected run only while the canvas still shows the graph that run was launched from
        if active and (active / "workflow.json").exists():
            try:
                if {n["type"] for n in json.loads((active / "workflow.json").read_text())["nodes"]} == {n["type"] for n in S["wf_nodes"]}: status = node_status(active)
            except Exception: pass
        val = _comp(palette=NODE_DEFS, saved=S["wf_nodes"], pvals=S["wf_pvals"], locked=False, status=status, preset=S["wf_preset"], height=max(430, vh - 128), key="wf_canvas", default=None)
        if val and val.get("nodes") is not None: S["wf_nodes"] = val["nodes"]; S["wf_pvals"] = val.get("pvals") or {}
    nodes = S["wf_nodes"] or []
    types = {n["type"] for n in nodes}

    with top[1]:
        st.markdown('<div class="ps-kick" style="margin:0 0 6px;line-height:1.4;font-size:1.05rem;letter-spacing:.08em">Inputs &amp; launch</div>', unsafe_allow_html=True)
        inp = st.container(key="wfin", height=max(340, vh - 170), border=False)
    with inp:
        exp_name = st.text_input("Experiment name", key="wf_exp", placeholder="e.g. WRN apo · covalent test", help="Names this run. It is kept with the results, so runs can be told apart later.")
        prot = st.file_uploader("Protein · PDB / mmCIF", type=["pdb", "cif", "mmcif"], help="Protein structure (PDB or mmCIF). Waters, ions and ligands are removed for docking.")
        lib = st.file_uploader("Compounds · CSV (id, smiles)", type=["smi", "txt", "csv", "tsv", "sdf", "xlsx", "xls"], help=FORMAT_HELP)
        if not lib: st.caption("Format: a CSV with a header row `id,smiles`, one compound per row, e.g. `aspirin,CC(=O)Oc1ccccc1C(=O)O`.")
        st.download_button("⬇ Download compound template (CSV)", TEMPLATE_CSV, file_name="pscreen_compounds_template.csv", mime="text/csv", key="tpl_dl",
                           help=FORMAT_HELP, use_container_width=True)
        kl = st.file_uploader("Known ligand · optional", type=["pdb", "sdf", "mol", "mol2"], disabled="known_ligand" not in types, help="Used by the Known-ligand site node: the site is defined from this co-crystal ligand.")
        exp = st.file_uploader("Experimental data · optional", type=["csv"], help="CSV with compound_id plus measured columns; joined to the results unchanged.")
        paths, info, nlib, receptor_pdb = {}, None, None, None
        tmp = RUNS / "_preview"; tmp.mkdir(parents=True, exist_ok=True)
        if prot:
            prot_path = tmp / hosting.safe_name(prot.name); prot_path.write_bytes(prot.getvalue())
            try: info = prepare_receptor(prot_path, tmp / "rec.pdb"); receptor_pdb = str(tmp / "rec.pdb")
            except Exception as e: st.error(f"Protein problem: {e}")
        lib_rep = None
        if lib:
            lp = tmp / hosting.safe_name(lib.name); lp.write_bytes(lib.getvalue())
            lib_rep = validate_library(lp)
            for e in lib_rep.errors: st.markdown(f'<div class="ps-err">✕ {e}</div>', unsafe_allow_html=True)
            if lib_rep.n_total and not lib_rep.errors:
                nlib = lib_rep.n_ok
                std = tmp / "library_standard.csv"; std.write_text(lib_rep.standard_csv()); paths["library_preview"] = str(std)
                st.caption(f"✓ {nlib} valid compounds of {lib_rep.n_total} rows" + (f" · {lib_rep.n_total - nlib} skipped (reported, never silently dropped)" if lib_rep.n_total > nlib else ""))
                with st.expander("Check how your compounds were read", expanded=bool(lib_rep.warnings)):
                    st.dataframe(lib_rep.table, hide_index=True, use_container_width=True, height=min(260, 38 + 35 * len(lib_rep.table)))
                if nlib > 1500 and not hosting.HOSTED: st.error("Library is larger than the supported ~1000-compound scale.")
        if prot and info: st.caption(f"✓ {info['n_residues']} residues · {info['n_chains']} chain(s)")
        tn = next((n for n in nodes if n["type"] == "target"), None)
        if tn and receptor_pdb:
            from ..residues import validate_target
            for r in parse_residues(tn["params"].get("residues", "")):
                v = validate_target(receptor_pdb, r)
                (st.caption(f"✓ {v['spec']} · nucleophile {v['nucleophile_atom']}") if v["ok"] else st.error(" ".join(v["errors"])))
        key = slug(exp_name)
        stage_dir = RUNS / "_staging"; stage_dir.mkdir(parents=True, exist_ok=True)
        ready_paths = {}
        for k2, f, pv in (("receptor", prot, None), ("library", lib, paths.get("library_preview")), ("known_ligand", kl, None), ("experimental_csv", exp, None)):
            if f:
                p = stage_dir / (Path(pv).name if pv else hosting.safe_name(f.name))
                p.write_bytes(Path(pv).read_bytes() if pv else f.getvalue()); ready_paths[k2] = str(p.resolve())
        cfg_kw, errs, warns = to_config(nodes, key, ready_paths, receptor_pdb)
        if not exp_name.strip(): errs = ["Name the experiment (field above).", *[e for e in errs if not e.startswith("Give the run")]]; cfg_kw = None
        for e in errs: st.markdown(f'<div class="ps-err">✕ {e}</div>', unsafe_allow_html=True)
        for w in warns: st.markdown(f'<div class="ps-warn">! {w}</div>', unsafe_allow_html=True)
        if cfg_kw and nlib:
            npk = max(1, len(cfg_kw["predictors"]) * cfg_kw["max_pockets_per_tool"] + bool(cfg_kw["known_ligand"]) + bool(cfg_kw["box_center"]))
            st.caption(f"≈ up to {nlib * npk} non-covalent docking jobs")
    with top[1]:
            if live_runs and cfg_kw:
                st.caption(f"{len(live_runs)} run(s) already active — runs are independent and share this machine's CPUs (each uses CPUs/job × Parallel).")
            cap_err = hosting.check_config(cfg_kw or {}, nlib, len(prot.getvalue()) if prot else 0, len(lib.getvalue()) if lib else 0,
                                           len(live_runs), len([r for r in list_runs() if not r.name.startswith("example")])) if cfg_kw else []
            for e in cap_err: st.markdown(f'<div class="ps-err">✕ {e}</div>', unsafe_allow_html=True)
            if hosting.HOSTED and not hosting.UNLIMITED and cfg_kw: st.caption(f"Free online limits: ≤{hosting.LIMITS['max_compounds']} compounds, ≤{hosting.LIMITS['max_residues']} covalent residues, 1 run at a time. Runs may wait in a queue.")
            if st.button("▶ Run analysis", type="primary", disabled=cfg_kw is None or bool(cap_err), use_container_width=True, key="wf_run"):
                cfg_kw = hosting.clamp_config(cfg_kw)
                rd = RUNS / key; (rd / "uploads").mkdir(parents=True, exist_ok=True)
                final = {}
                for k2 in ("receptor", "library", "known_ligand", "experimental_csv"):
                    if cfg_kw.get(k2):
                        dst = rd / "uploads" / Path(cfg_kw[k2]).name; dst.write_bytes(Path(cfg_kw[k2]).read_bytes()); final[k2] = str(dst.resolve())
                cfg_kw.update(final, name=exp_name.strip()); RunConfig(**cfg_kw).save(rd)
                (rd / "workflow.json").write_text(json.dumps({"nodes": nodes, "pvals": S["wf_pvals"]}, indent=1))
                launch(rd); S["active_run"] = key; go("Runs")
