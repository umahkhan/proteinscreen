"""pscreen · Pockets — find candidate binding pockets (conventional and cryptic) on a protein structure. Pocket finding only: no docking.

Light enough for a free Streamlit host: P2Rank (Java), Lacuna (pip) and, where installed, PocketMiner. Every visitor is a separate Streamlit session with
its own temporary folder; nothing is stored after the session (old folders are deleted).   Run locally:  streamlit run pockets_app.py
"""
import os
import shutil
import tempfile
import time
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components

os.environ.setdefault("PSCREEN_AUTO_TOOLS", "1")                 # fetch P2Rank on first use when it is not installed
from pscreen import hosting, pocket_tools, pocket_report as pr
from pscreen.pipeline import prepare_receptor
from pscreen.ui import style, viewer3d

st.set_page_config(page_title="pscreen · Pockets", page_icon="🧬", layout="wide")
style.inject()
BASE = Path(tempfile.gettempdir()) / "pscreen_pockets"; BASE.mkdir(exist_ok=True)
for old in BASE.iterdir():                                         # clean-up: sessions older than 2 h
    try:
        if time.time() - old.stat().st_mtime > 7200: shutil.rmtree(old, ignore_errors=True)
    except OSError: pass
S = st.session_state
if "work" not in S: S["work"] = tempfile.mkdtemp(dir=BASE)
WORK = Path(S["work"])


def chain_ids(path: Path) -> list[str]:
    if path.suffix.lower() not in (".pdb", ".ent"): return []
    return sorted({l[21] for l in path.read_text(errors="ignore").splitlines() if l.startswith("ATOM") and len(l) > 21 and l[21].strip()})


def keep_chains(src: Path, chains: list[str]) -> Path:
    out = WORK / f"chains_{''.join(chains)}.pdb"
    out.write_text("".join(l for l in src.read_text(errors="ignore").splitlines(True) if not l.startswith(("ATOM", "HETATM", "ANISOU", "TER")) or (len(l) > 21 and l[21] in chains)))
    return out


st.markdown(f'<div class="ps-brand">{style.LOGO}<div><div class="ps-brand-t">pscreen · Pockets</div><div class="ps-brand-s">Find candidate binding pockets, including cryptic ones</div></div></div>', unsafe_allow_html=True)
st.markdown('<div class="ps-div" style="margin:.5rem 0 .8rem"></div>', unsafe_allow_html=True)
st.caption("Upload a protein structure, choose the pocket models, and explore the predicted pockets in 3D. Free, no account. Your file is used only for this session and is deleted afterwards. "
           "Predictions only: a predicted pocket is not evidence that a compound binds or that the site is druggable.")

left, right = st.columns([1, 2.4], gap="large")
with left:
    st.markdown('<div class="ps-kick">1 · Protein structure</div>', unsafe_allow_html=True)
    how = st.radio("Source", ["Upload a file", "Fetch from the PDB"], horizontal=True, label_visibility="collapsed")
    src = None
    EX = Path(__file__).parent / "examples_data" / "9S1A_chainA_apo.pdb"
    if how == "Upload a file":
        up = st.file_uploader("PDB or mmCIF file", type=["pdb", "ent", "cif", "mmcif"], help="A protein structure (PDB or mmCIF). Waters, ions and ligands are removed.")
        if up:
            S["example"] = False
            if len(up.getvalue()) > pr.LIMITS["max_mb"] * 2**20: st.error(f"File is larger than {pr.LIMITS['max_mb']} MB.")
            else: src = WORK / hosting.safe_name(up.name); src.write_bytes(up.getvalue())
        elif EX.exists():
            if st.button("Try the WRN helicase example", use_container_width=True, help="Apo WRN helicase (PDB 9S1A, chain A), the target of covalent drugs at Cys727. A real structure, so you can see typical results."):
                S["example"] = True
            if S.get("example"):
                src = WORK / EX.name; shutil.copy(EX, src); st.caption(f"Using the example structure: {EX.name}")
    else:
        pid = st.text_input("PDB ID", placeholder="e.g. 3MXF", max_chars=4)
        if pid:
            try: src = pr.fetch_pdb(pid, WORK)
            except Exception as e: st.error(str(e))
    chains = chain_ids(src) if src else []
    sel_chains = st.multiselect("Chains", chains, default=chains[:1], help="Choose the chain(s) to analyse. Using one chain is faster and avoids duplicate pockets.") if len(chains) > 1 else chains

    st.markdown('<div class="ps-kick" style="margin-top:.8rem">2 · Pocket models</div>', unsafe_allow_html=True)
    use_p2 = st.checkbox("P2Rank — conventional pockets", True, help="Fast geometry + machine-learning pocket predictor (seconds).")
    use_lac = st.checkbox("Lacuna — cryptic pockets (conformational ensemble)", True, disabled=pocket_tools.lac_python() is None,
                          help="Generates alternative conformations and looks for pockets that open in them. Takes a few minutes.")
    n_conf = st.slider("Lacuna conformers", 5, 20, 10, disabled=not use_lac, help="More conformers explore more shapes but take longer.")
    pm_ok = pocket_tools.pm_possible()
    use_pm = st.checkbox("PocketMiner — cryptic pocket propensity", False, disabled=not pm_ok,
                         help=("Predicts which residues are likely to form cryptic pockets (original PocketMiner model). The first use on a server sets up TensorFlow and takes a few minutes."
                               if pm_ok else "Not available on this server; available in the local install of pscreen."))
    run = st.button("▶ Find pockets", type="primary", use_container_width=True, disabled=src is None or not (use_p2 or use_lac or use_pm))

if run and src:
    S.pop("pk", None)
    with right, st.status("Finding pockets…", expanded=True) as status_box:
        try:
            use = keep_chains(src, sel_chains) if sel_chains and chains else src
            st.write("Preparing the structure…")
            rec = WORK / "receptor.pdb"; info = prepare_receptor(use, rec)
            if info["n_residues"] > pr.LIMITS["max_residues"]: raise ValueError(f"{info['n_residues']} residues is above the limit of {pr.LIMITS['max_residues']}. Choose fewer chains.")
            if use_lac and info["n_residues"] > pr.LIMITS["lacuna_max_residues"]:
                use_lac = False; st.warning(f"Lacuna was skipped: {info['n_residues']} residues is above its limit of {pr.LIMITS['lacuna_max_residues']} on this server. Choose one chain.")
            pockets, stat = [], {}
            for tool, on, fn in (("p2rank", use_p2, lambda: pocket_tools.run_p2rank(str(rec), WORK / "p2rank", 5, padding=4.0)),
                                 ("lacuna", use_lac, lambda: pocket_tools.run_lacuna(str(rec), WORK / "lacuna", n_conformers=n_conf, top=5)),
                                 ("pocketminer", use_pm, lambda: pocket_tools.run_pocketminer(str(rec), WORK / "pocketminer", max_pockets=5)[0])):
                if not on: continue
                st.write(f"Running {pr.MODEL_OF[tool]}…"); t0 = time.time()
                try:
                    if tool == "pocketminer": pocket_tools.ensure_pocketminer_env(st.write)
                    pk = fn(); pockets += pk; stat[tool] = f"{len(pk)} pocket(s) in {time.time() - t0:.0f} s" if pk else "ran, no pockets found"
                except Exception as e: stat[tool] = f"FAILED: {e}"; st.warning(f"{pr.MODEL_OF[tool]} failed: {str(e)[:300]}")   # never silently replaced by another tool
            if not pockets: raise RuntimeError("No pockets were found by the selected models.")
            panel = viewer3d.panel(src.stem, f"{len(pockets)} predicted site(s) on the input structure", str(rec), [pr.pocket_dict(p) for p in pockets], [], None, None, pockets[0].pocket_id)
            S["pk"] = {"name": src.stem, "rec": str(rec), "pockets": pockets, "status": stat, "panel": panel, "info": info}
            status_box.update(label="Done", state="complete", expanded=False)
        except Exception as e:
            status_box.update(label="Could not finish", state="error"); st.error(str(e))

with right:
    R = S.get("pk")
    if not R:
        st.markdown('<div class="ps-empty"><b>No pockets yet</b><br>Add a structure and press <i>Find pockets</i>. Try PDB ID <code>3MXF</code> (BRD4) for a quick example.</div>', unsafe_allow_html=True)
    else:
        st.markdown(f'<div class="ps-kick">Results · {R["name"]} · {R["info"]["n_residues"]} residues</div>', unsafe_allow_html=True)
        for t, v in R["status"].items(): st.caption(f"{pr.MODEL_OF[t]}: {v}")
        components.html(viewer3d.explorer_html([R["panel"]]), height=620)
        tbl = pr.pockets_table(R["pockets"])
        st.dataframe(tbl, hide_index=True, use_container_width=True, column_config={"lining residues": st.column_config.TextColumn(width="large")})
        st.download_button("⬇ Download results (zip: CSV, JSON, structure, PyMOL scene)", pr.bundle(R["rec"], R["pockets"], R["status"], R["panel"]),
                           file_name=f"pscreen_pockets_{R['name']}.zip", mime="application/zip", type="primary")
        st.caption("Pockets are computational predictions, not experimentally confirmed. Cryptic-pocket predictions on closed (apo) structures are hypotheses. "
                   "Tools: P2Rank, Lacuna, PocketMiner. Please cite them.")

with st.expander("About pscreen · Pockets"):
    st.markdown("""
**What it does.** Predicts candidate binding pockets on a protein structure: conventional pockets with **P2Rank**, cryptic pockets with **Lacuna** (an ensemble of
conformations) and **PocketMiner** (per-residue cryptic-pocket probability). It does not dock compounds.

**How to read the results.** Pockets are computational predictions. They are not experimentally confirmed, and they do not show that a compound binds or that the site is druggable.
Predictions on closed (apo) structures are hypotheses. Tool scores from different models are shown side by side and never combined.

**Privacy.** Your file is used only in your session and deleted afterwards.

**Cite the tools.** P2Rank (Krivák & Hoksza, 2018) · PocketMiner (Meller et al., Nature Communications 2023) · Lacuna.
""")
