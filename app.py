"""pscreen — structure-first local app: define a screen, monitor it, explore sites and compounds in 3D beside a compact ranked table, export."""
import io
import os
import subprocess
import zipfile
import sys
from pathlib import Path
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem

from pscreen.docking import MODE_COVALENT, MODE_NONCOVALENT
from pscreen.results import best_pose_per_compound
from pscreen.runio import load_run
from pscreen.validation.retro import evaluate
from pscreen.viewer import pose_block
from pscreen import pymol_export
from pscreen.ui import style, viewer3d
from pscreen.ui.components import SIGNAL_LABEL, build_scene, histogram, size_scatter, top_bar, warhead_bar
from pscreen.ui.dialogs import compound_dialog
from pscreen.ui.tables import ADMET_COMPACT, EXTRA_COLUMNS, add_chem_columns, build_table, has_wh, structure_grid
from pscreen.ui import target_picker, workflow
from pscreen.ui.layout import report_viewport, viewport_h
from pscreen import hosting

st.set_page_config(page_title="pscreen", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")
style.inject()
if hosting.HOSTED:                                  # private per-visitor workspace, addressed by ?ws=<token> (bookmark it to come back)
    _ws = st.query_params.get("ws")
    if not hosting.valid_token(_ws): _ws = hosting.new_token(); st.query_params["ws"] = _ws
    st.session_state["ws"] = _ws; hosting.ensure_workspace(_ws); hosting.touch(_ws)
    @st.cache_resource
    def _janitor(): return hosting.purge_old()
    _janitor()
RUNS = workflow.runs_root()
SIGNALS = ["gnina_cnn_vs", "gnina_cnn_pose_score", "gnina_cnn_affinity_pK", "gnina_minimized_affinity_kcal"]


def pid_alive(pid: int) -> bool:
    try: os.kill(pid, 0); return True
    except OSError: return False


def launch(rundir: Path):
    log = open(rundir / "stdout.log", "ab")
    p = subprocess.Popen([sys.executable, "-m", "pscreen.pipeline", str(rundir)], stdout=log, stderr=log, start_new_session=True)
    (rundir / "pid").write_text(str(p.pid))


def conformer_text(ensemble_pdb: str, idx: int) -> str:
    out, on = [], False
    for l in open(ensemble_pdb):
        if l.startswith("MODEL"): on = int(l.split()[1]) - 1 == idx
        elif l.startswith("ENDMDL"): on = False
        elif on and l.startswith(("ATOM", "HETATM", "TER")): out.append(l)
    return "".join(out) + "END\n"


report_viewport()
_g = st.session_state.pop("_goto", None)
if "_deeplink_done" not in st.session_state:        # open straight to a page / experiment:  http://localhost:8501/?page=Results&run=<run folder name>
    st.session_state["_deeplink_done"] = True
    _qp = st.query_params
    if _qp.get("page") in ("Workflow", "Runs", "Results"):
        _g = {**(_g or {}), "nav": _qp["page"], **({"res_run": _qp["run"]} if _qp.get("run") and (workflow.runs_root() / _qp["run"]).exists() else {})}
st.session_state.setdefault("nav", "Workflow")
if _g: st.session_state.update(_g)                  # navigation requested by a button on the previous run (applied before the nav widget exists)
_n_live = sum(workflow.run_state(r)[0] in workflow.ACTIVE for r in workflow.list_runs())
RANK_HELP = ("**Pick one score to sort by; scores are never combined.**\n\n"
             "- **CNN VS** – pose confidence × predicted affinity. Good default: favours compounds that are both well placed and predicted to bind well.\n"
             "- **CNN pose** – how realistic the docked pose looks (0–1, higher is better).\n"
             "- **CNN affinity** – the network's estimate of binding strength (pKd-like, higher is better).\n"
             "- **Vina-like** – classic docking energy in kcal/mol (more negative is better); blind to pose quality.\n\n"
             "All are model predictions for prioritising compounds, not evidence of binding or potency.")
PAGES = {"Workflow": "1 · Workflow", "Runs": "2 · Runs", "Results": "3 · Results"}
nh1, nh2, nh3 = st.columns([1.5, 2.6, 2.2], vertical_alignment="center")
nh1.markdown(f'<div class="ps-brand">{style.LOGO}<div><div class="ps-brand-t">pscreen</div><div class="ps-brand-s">Pocket &amp; compound prioritization</div></div></div>', unsafe_allow_html=True)
page = nh2.segmented_control("Page", list(PAGES), key="nav", format_func=PAGES.get, label_visibility="collapsed") or "Workflow"
_live_txt = (f'<span class="dot"></span>{_n_live} experiment' + ("s" if _n_live != 1 else "") + " running") if _n_live else "No experiments running"
nh3.markdown(f'<div class="ps-live">{_live_txt}</div>', unsafe_allow_html=True)
st.markdown('<div class="ps-div" style="margin:.35rem 0 .7rem"></div>', unsafe_allow_html=True)
if hosting.HOSTED:
    _base = os.environ.get("PSCREEN_PUBLIC_URL", "").rstrip("/")
    with st.expander("🔒 Your private workspace: bookmark this link to come back to your experiments", expanded=False):
        st.code(f"{_base}/?ws={st.session_state['ws']}" if _base else f"?ws={st.session_state['ws']}  (add this to the page address)", language=None)
        st.caption(f"Other visitors cannot see your runs, but your files are processed and stored on the computer of whoever operates this server until they are deleted ({int(hosting.RETENTION_H)} hours without use). "
                   "Do not upload confidential structures or compounds unless you trust the operator. No account and no tracking. "
                   "Predictions only: docking scores and pockets are not evidence of binding, engagement, potency or selectivity.")

if page == "Workflow":
    workflow.render()
    st.stop()
if page == "Runs":
    workflow.runs_page()
    st.stop()

# =============================================================================================================== RESULTS
else:
    runs = [r for r in workflow.list_runs() if (r / "results.sqlite").exists()]
    if not runs:
        st.markdown('<div class="ps-empty"><b>No finished experiments yet</b><br>Results appear here once a run in <i>2 · Runs</i> has finished.</div>', unsafe_allow_html=True); st.stop()
    _names = [r.name for r in runs]
    if st.session_state.get("res_run") not in _names: st.session_state["res_run"] = _names[0]
    sw0, sw1, _sw = st.columns([0.55, 1.6, 2.6], vertical_alignment="center", gap="small")
    sw0.markdown('<div class="ps-kick" style="white-space:nowrap">Switch experiment</div>', unsafe_allow_html=True)
    sw1.selectbox("Switch experiment", _names, key="res_run", label_visibility="collapsed", format_func=lambda n: workflow.run_title(workflow.runs_root() / n), help="Everything below belongs to the experiment selected here.")
    st.markdown('<div class="ps-div" style="margin:.3rem 0 .9rem"></div>', unsafe_allow_html=True)
    rd = workflow.runs_root() / st.session_state["res_run"]
    if st.session_state.get("_res_loaded") != rd.name:                      # never carry one experiment's selections into another
        for _k in ("pm_models", "conf_sel", "dbl_row", "_dbl_done", "xc", "section", "rev_target_mode", "rev_target_type", "rev_target_res", "rev_target_txt"): st.session_state.pop(_k, None)
        st.session_state["_res_loaded"] = rd.name
    _st, _d = workflow.run_state(rd)
    _ctx_pending = True
    st.markdown(f'<div class="ps-ctx ps-ctx2"><div class="ps-ctx2-l"><span class="ps-kick">Viewing results for</span>'
                f'<div class="ps-big">{workflow.run_title(rd)} {workflow.chip(_st)}</div><div class="ps-nodes">{workflow.run_nodes_html(rd)}</div></div>'
                f'<div class="ps-ctx2-r"><div class="ps-fact"><span>Run ID</span><code>{rd.name}</code></div>'
                f'<div class="ps-fact"><span>Started</span><b>{workflow.started_str(rd, _d)}</b></div></div></div>', unsafe_allow_html=True)
    if _st != "finished":                                                   # partial runs have no complete result set yet
        _msg = {"running": "This experiment is still running. Results appear when it finishes (follow progress in 2 · Runs).", "queued": "This experiment is waiting in the queue.",
                "failed": "This experiment failed, so there are no results. See the log in 2 · Runs.", "interrupted": "This experiment was interrupted. Resume it in 2 · Runs."}.get(_st, "No results yet.")
        st.info(_msg); st.stop()

    @st.cache_resource(show_spinner="Loading screen…")
    def cached(rd_s, mtime, target):
        return load_run(Path(rd_s), target)
    R0 = cached(str(rd), (rd / "results.sqlite").stat().st_mtime, None)
    with st.sidebar:
        st.markdown("**Target residue**")
        first_c = next((p.center for p in R0["pockets"]), None)
        tgt_spec, _ = target_picker.pick(R0["receptor"], "rev_target", R0["docked_target"], first_c)
        if tgt_spec != R0["docked_target"]:
            st.caption("Geometry columns for the noncovalent lane are recomputed for this residue. The covalent lane stays as docked" + (f" ({R0['docked_target']})." if R0["docked_target"] else " (none)."))
        st.caption(f"GNINA `{str(R0['meta'].get('gnina_version'))[:12]}` · seed {R0['meta'].get('config', {}).get('seed')}")
    R = R0 if tgt_spec == R0["docked_target"] else cached(str(rd), (rd / "results.sqlite").stat().st_mtime, tgt_spec)
    res, pockets, prepared, meta, target = R["results"], R["pockets"], R["prepared"], R["meta"], tgt_spec
    docked = res[res.processing_status == "docked"]


    adv = st.sidebar.checkbox("Show advanced sections", False, help="Overview, Pockets & conformations, Experimental data, QC, Export (hidden by default).")
    SECTIONS = ["Explore"] + (["Overview", "Pockets & conformations", "Experimental data", "Run status & QC", "Export"] if adv else [])
    section = (st.segmented_control("Section", SECTIONS, default="Explore", key="section", label_visibility="collapsed") or "Explore") if adv else "Explore"

    ACTIVE_IDS = [p.pocket_id for p in pockets]

    def lens():
        c1, c2, c3, c4, c5 = st.columns([1.1, 1.5, 0.55, 1.2, 0.55])
        has_cov = (docked.docking_mode == MODE_COVALENT).any()
        site = c1.selectbox("Rank compounds at", ["All sites"] + [p.pocket_id for p in pockets], help="'All sites': each compound appears once, at its best site. All sites stay mapped in the explorer either way.")
        sig = c2.selectbox("Rank by", SIGNALS, label_visibility="collapsed", format_func=lambda s: SIGNAL_LABEL[s], help=RANK_HELP)
        topn = c3.number_input("Top N", 3, 200, 15)
        lane = c4.segmented_control("Docking lane", ([MODE_COVALENT] if has_cov else []) + [MODE_NONCOVALENT], default=MODE_COVALENT if has_cov else MODE_NONCOVALENT,
                                    format_func=lambda m: "Noncovalent" if m == MODE_NONCOVALENT else "Covalent (adduct)", help="Scores from different lanes are never compared.") or MODE_NONCOVALENT
        with c5.popover("⚙", help="Layout"):
            layout = st.radio("Layout", ["Explorer above table", "Explorer beside table", "Explorer only", "Table only"])
            height = st.slider("Explorer height (px)", 700, 1400, 940, 20)
        return lane, site, sig, int(topn), layout, height

    def rank_ids():
        """Compounds are only docked at some sites; rank at the selected sites that were docked (else at all docked sites, e.g. when only a display-only model is selected)."""
        dk = set(res.pocket_id.dropna()); ids = [i for i in ACTIVE_IDS if i in dk]
        return ids or [i for i in dk]

    def ranked(lane, site, sig, topn=None):
        """ONE row per compound: its best pose over the selected site(s) by the chosen signal. Scores at the other sites are kept in `site_scores`."""
        asc = sig == "gnina_minimized_affinity_kcal"
        b = best_pose_per_compound(res, sig, lane)
        b = b[b.pocket_id.isin(rank_ids())]
        if site != "All sites": b = b[b.pocket_id == site]
        if b.empty: return pd.DataFrame()
        b = b.sort_values(sig, ascending=asc)
        ss = b.groupby("ligand_key").apply(lambda g: " · ".join(f"{r.pocket_id} {r[sig]:.2f}" for _, r in g.iterrows()), include_groups=False).to_dict()
        out = b.drop_duplicates("ligand_key", keep="first").copy()
        out["rank"] = range(1, len(out) + 1); out["top10"] = out["rank"] <= max(1, int(round(0.10 * len(out)))); out["site_scores"] = out.ligand_key.map(ss)
        return out

    def pocket_dict(p): return {"id": p.pocket_id, "residues": p.residues, "center": list(p.center), "size": list(p.size), "source": p.source, "residue_scores": p.metrics.get("residue_scores"), "residue_score_label": p.metrics.get("residue_score_label", "")}

    def lac_vol(q):
        return {int(k): float(v) for k, v in (q.metrics.get("volume_by_conformer") or {}).items() if v}

    def conformations(lacs):
        """ALL Lacuna conformers (0 = input). Each Lacuna site has its own per-conformer volumes; a conformer lists the sites detected in it (ranked by largest volume)."""
        vbs = {q.pocket_id: lac_vol(q) for q in lacs}
        n = max(int(q.metrics.get("n_conformers") or 0) for q in lacs) or 1
        docked = {q.pocket_id: q.metrics.get("open_conformer_idx") for q in lacs}
        opts = {"Input structure (conformer 0)": ("input", None)}
        have = sorted({i for v in vbs.values() for i in v if i}, key=lambda i: -max(v.get(i, 0) for v in vbs.values()))
        for rank, idx in enumerate(have, 1):
            here = [(pid, v[idx]) for pid, v in vbs.items() if idx in v]
            tag = " · ".join(f"{pid} {vol:.0f} Å³" for pid, vol in here)
            dk = [pid for pid, _ in here if docked.get(pid) == idx]
            opts[f"#{rank} · conformer {idx} · {tag}" + (f" · docked: {', '.join(dk)}" if dk else "")] = ("conf", idx)
        for idx in range(1, n):
            if idx not in have: opts[f"conformer {idx} · no pocket detected"] = ("conf", idx)
        return opts, docked

    # ====================================================================== EXPLORE
    if section == "Explore":
        def model_of(p):
            t = str(p.metrics.get("tool", ""))
            return "Lacuna" if t.startswith("Lacuna") else "P2Rank" if t.startswith("P2Rank") else "PocketMiner" if t.startswith("PocketMiner") else \
                   "Known ligand" if p.source == "known_ligand" else "User box" if p.source == "user_box" else p.source.replace("_", " ").title()
        models = list(dict.fromkeys(model_of(p) for p in pockets)); has_cov = (docked.docking_mode == MODE_COVALENT).any(); rec0 = R["receptor"]
        explorer_box, ctl_box = st.columns([1.15, 1], gap="medium")
        with explorer_box:
            h0, h1 = st.columns([0.55, 1], vertical_alignment="center")
            h0.markdown(f'<div class="ps-hd" style="margin:0"><span class="ps-kick" style="white-space:nowrap">3D explorer</span><span class="ps-s" style="white-space:nowrap">Target {R["docked_target"] or "—"}</span></div>', unsafe_allow_html=True)
            with h1:
                chosen = st.pills("Pocket model", models, selection_mode="multi", default=[models[0]], key="pm_models", label_visibility="collapsed", help="Which pocket-prediction model(s) to map on the protein. Compounds are ranked at the sites of the selected model(s) only.")
                chosen = chosen or [models[0]]
            pockets_sel = [p for p in pockets if model_of(p) in chosen]; ACTIVE_IDS[:] = [p.pocket_id for p in pockets_sel]
            lac = [p for p in pockets_sel if model_of(p) == "Lacuna"]; lp = lac[0] if lac else None; conf = ("input", None); compare = False
            if lp is not None:
                opts, _ = conformations(lac); labels = list(opts); default_label = next((l for l in labels if "docked" in l), labels[0])
                conf = opts[st.selectbox("Conformation (Lacuna)", labels, index=labels.index(default_label), key="conf_sel", label_visibility="collapsed", help="Lacuna conformations ranked by pocket volume; the best-ranked one was used for docking. Everything is drawn on the one you pick.")]
        with ctl_box:
            t0, t1 = st.columns([5, 1], vertical_alignment="center")
            t0.markdown('<div class="ps-h" style="margin:0">Ranked compounds</div>', unsafe_allow_html=True)
            with t1.popover("⚙", help="Layout"):
                layout = st.radio("Layout", ["Side by side", "Explorer only", "Table only"]); ex_height = st.slider("Table height (px)", 200, 1000, int(min(1000, max(180, viewport_h() - 470))), 10)
                compare = st.checkbox("Compare with input structure", False, disabled=lp is None, help="Two camera-linked panels: closed input structure vs the chosen conformation.")
            ctl = st.container(key="ctl")
            c4, c1, c2, c3 = ctl.columns([1.9, 1.05, 1.6, 0.85], gap="small")
            site = c1.selectbox("Site", options=["All selected sites"] + rank_ids(), format_func=lambda x: "All" if x.startswith("All") else x, help="'All selected sites': each compound appears once, at its best site.")
            site = "All sites" if site.startswith("All") else site
            sig = c2.selectbox("Rank by", SIGNALS, format_func=lambda s_: SIGNAL_LABEL[s_].split(" (")[0].replace(" score","").replace("-like score","-like"), help=RANK_HELP)
            topn = int(c3.number_input("Top N", 3, 200, 15, help="Number of top compounds to list"))
            lane = c4.segmented_control("Docking lane", ([MODE_COVALENT] if has_cov else []) + [MODE_NONCOVALENT], default=MODE_COVALENT if has_cov else MODE_NONCOVALENT,
                                        format_func=lambda m: "Non-cov." if m == MODE_NONCOVALENT else "Covalent", help="Docking lane. Covalent: the compound is docked in its bonded (adduct) form attached to the target residue; only poses with the bond satisfied are ranked. Non-covalent: the intact compound is docked freely. Scores from the two lanes are never compared.") or (MODE_COVALENT if has_cov else MODE_NONCOVALENT)
        best = ranked(lane, site, sig)
        if best.empty: st.info("No docked poses for the selected pocket model(s) in this lane."); st.stop()
        tab = build_table(R, best, topn, len(pockets_sel) > 1)
        if "Site" in tab: tab = tab.rename(columns={"Site": "Best site"})
        key = f"tbl_{lane}_{site}_{sig}_{topn}_{'-'.join(rank_ids())}_{rd.name}"
        prev = (st.session_state.get(key) or {}).get("selection", {}).get("rows", []) if isinstance(st.session_state.get(key), dict) else []
        primary_key = tab.iloc[prev[0]]["_ligand_key"] if prev and prev[0] < len(tab) else tab.iloc[0]["_ligand_key"]
        admet_cols = [c for _, c, _ in ADMET_COMPACT if c in tab.columns]

        PML = {"panel": None, "slot": None}
        def _fill_pml():
            if PML["slot"] is None or PML["panel"] is None: return
            PML["slot"].download_button("⬇ PyMOL (.pml)", pymol_export.bundle(PML["panel"]), f"pscreen_pymol_{rd.name}.zip", mime="application/zip", key=f"pml_{rd.name}",
                                        help="Zip with a PyMOL script in the explorer's design (pockets, compounds, covalent target, scenes F1–F4) plus the structure files it loads. Predictions only.")
        def show_table(height):
            view_mode = "Table"
            more = st.multiselect("More columns", admet_cols + EXTRA_COLUMNS, key="xc", placeholder="More columns: ADMET (solubility, hERG, AMES, DILI), MW, cLogP, TPSA, Lipinski, scaffold, similarity …", label_visibility="collapsed") if view_mode == "Table" else []
            tt = add_chem_columns(R, tab, [c for c in more if c in EXTRA_COLUMNS], primary_key)
            cfg = {"Rank": st.column_config.NumberColumn(width="small"), "Compound": st.column_config.ImageColumn("Compound", width="medium", help="Warhead atoms shaded; reactive atom in red"),
                   "Warhead structure": st.column_config.ImageColumn("Warhead", width="small"), "Name": st.column_config.TextColumn("ID", width="medium"),
                   "CNN pose": st.column_config.ProgressColumn("CNN pose", min_value=0, max_value=1, format="%.2f", help="Probability the pose resembles a good crystal-like pose"),
                   "Vina-like": st.column_config.NumberColumn("Vina-like", format="%.2f", help="kcal/mol; more negative = better"), "CNN pK": st.column_config.NumberColumn(format="%.2f"),
                   "Reactive→target (Å)": st.column_config.NumberColumn("Reactive→target Å", format="%.1f"),
                   "Solubility logS": st.column_config.NumberColumn("Sol. logS", format="%.1f", help="ADMET-AI aqueous solubility, log(mol/L)"),
                   "hERG risk": st.column_config.ProgressColumn("hERG p", min_value=0, max_value=1, format="%.2f", help="ADMET-AI probability, not calibrated risk"),
                   "AMES": st.column_config.ProgressColumn("AMES p", min_value=0, max_value=1, format="%.2f"), "DILI": st.column_config.ProgressColumn("DILI p", min_value=0, max_value=1, format="%.2f"),
                   "Best site": st.column_config.TextColumn("Best site", help="Where this compound scores best by the selected signal."),
                   "Similarity to selected": st.column_config.ProgressColumn("Similarity to selected", min_value=0, max_value=1, format="%.2f", help="Morgan radius-2 Tanimoto vs the selected compound")}
            order = [c for c in ["Rank", "Compound", "Name", "Best site", "CNN pose", "Vina-like", "Warhead structure", "Warhead", "Reactive→target (Å)"] + more if c in tt.columns]
            if view_mode == "Grid":
                st.image(structure_grid(R, tt), use_container_width=True)
                pick = st.multiselect("Compounds to overlay in the explorer (up to 6)", list(tt["Name"]), default=[], max_selections=6, key=f"gridpick_{key}")
                rows_ = [int(tt.index[tt["Name"] == n][0]) for n in pick]
            else:
                ev = st.dataframe(tt[order + ["_pocket", "_ligand_key"]], hide_index=True, use_container_width=True, column_config=cfg, column_order=order, on_select="rerun", selection_mode="multi-row", key=key, height=height, row_height=70)
                rows_ = ev.selection.rows
            c1, cP, c2 = st.columns([0.8, 1.15, 2.4], vertical_alignment="center"); PML["slot"] = cP.empty()
            if PML["panel"] is not None: _fill_pml()
            exp_cols = ["rank", "compound_id", "canonical_smiles", "pocket_id", "site_scores", "gnina_cnn_pose_score", "gnina_minimized_affinity_kcal", "gnina_cnn_affinity_pK", "gnina_cnn_vs", "warhead_classes", "reactive_to_target_distance", "covalent_geometry", "mol_weight", "clogp", "tpsa"] + [c for c in best.columns if str(c).startswith(("chemprop__", "exp_"))]
            c1.download_button("⬇ CSV", best.head(topn)[[c for c in exp_cols if c in best.columns]].to_csv(index=False), "pscreen_ranked_compounds.csv", help="One row per compound: all signals and any ADMET / experimental columns (no images).")
            c2.caption("Tick rows to overlay · double-click a row for details · Docking scores are predictions, not evidence of binding, engagement or potency.")
            return rows_

        if layout == "Table only":
            with ctl_box: rows = show_table(ex_height + 120)
            colv = None
        elif layout == "Explorer only": colv = explorer_box; rows = None
        else:
            colv = explorer_box
            with ctl_box: rows = show_table(ex_height)
        holder = None
        sel_rows = prev if rows is None else rows
        sel = [tab.iloc[i] for i in sel_rows if i < len(tab)][:6]          # nothing is overlaid until the user selects a compound
        items = []
        for ci, x in enumerate(sel):          # every selected compound at EVERY selected site where it docked
            for p in [q for q in pockets if q.pocket_id in rank_ids()]:
                sc = build_scene(R, [x["Name"]], p, lane, target, {})
                if sc: it = sc[0]; it["site"] = p.pocket_id; it["color_idx"] = ci; it["label"] = x["Name"]; items.append(it)

        if colv is not None:
            with colv:
                def conf_receptor():
                    if lp is None: return rec0
                    if conf[0] == "input": return rec0
                    ens = lp.metrics.get("ensemble_pdb")
                    if conf[1] == lp.metrics.get("open_conformer_idx") and lp.receptor_path: return lp.receptor_path
                    if ens and Path(ens).exists():
                        cf = rd / "renders" / f"conf_{lp.pocket_id}_{conf[1]}.pdb"; cf.parent.mkdir(exist_ok=True); cf.write_text(conformer_text(ens, int(conf[1]))); return str(cf)
                    return rec0
                prec = conf_receptor()
                # a Lacuna site exists only in the conformers where it was detected: never draw it on a conformer where it was not
                cidx = 0 if conf[0] == "input" else int(conf[1])
                shown = [p for p in pockets_sel if model_of(p) != "Lacuna" or cidx in lac_vol(p)]
                all_p = [pocket_dict(p) for p in shown]
                refp = next(((p.metrics or {}).get("reference_ligand") for p in pockets_sel if (p.metrics or {}).get("reference_ligand")), None)
                vtxt = lambda v: "pocket not detected here" if not v else f"pocket volume {v:.0f} Å³"
                title = "Input structure" if prec == rec0 else (f"Lacuna conformation {conf[1]}" if conf[1] is not None else "Lacuna conformation")
                subtitle = " · ".join(chosen) + (f" · {len(shown)} site(s)" if shown else " · no pocket detected in this conformer") + (f" · {len(sel)} compound(s)" if sel else " · select a compound to overlay it")
                ens = lp.metrics.get("ensemble_pdb") if lp is not None else None
                einfo = {"n": int(lp.metrics.get("n_conformers") or 0), "open": [int(k) for k in (lp.metrics.get("volume_by_conformer") or {})], "docked": lp.metrics.get("open_conformer_idx")} if ens else None
                B = viewer3d.panel(title, subtitle, prec, all_p, items, target, refp, pockets_sel[0].pocket_id if pockets_sel else None, ensemble_path=ens, ens_info=einfo)
                panels = [B]
                PML["panel"] = B; _fill_pml()
                if compare and lp is not None:
                    A = viewer3d.panel("Input structure (closed)", f"{lp.pocket_id} · {vtxt(lp.metrics.get('volume_input_A3'))}", rec0, all_p, [], target, None, lp.pocket_id); panels = [A, B]
                components.html(viewer3d.explorer_html(panels), height=700)
        if holder is not None:
            with holder: show_table(560)
        # double-click a row -> details: a script maps the click to a row index and writes "row:nonce" into a hidden text box
        st.markdown("<style>.st-key-dbl_row{position:absolute;width:0;height:0;overflow:hidden;opacity:0}</style>", unsafe_allow_html=True)
        dv = st.text_input("dbl", key="dbl_row", label_visibility="collapsed")
        if dv and dv != st.session_state.get("_dbl_done"):
            st.session_state["_dbl_done"] = dv
            try: ri = int(dv.split(":")[0])
            except ValueError: ri = -1
            if 0 <= ri < len(tab):
                x_ = tab.iloc[ri]
                pk = next(p for p in pockets if p.pocket_id == x_["_pocket"]); r_ = best[best.compound_id == x_["Name"]].iloc[0]
                sc = build_scene(R, [x_["Name"]], pk, lane, target, {})
                if sc: compound_dialog(R, sc[0], lane, bool(r_.top10), target)
        components.html("""<script>
          const P = window.parent; if (!P.__psDbl) { P.__psDbl = true;
            P.document.addEventListener('dblclick', e => { const df = e.target.closest('[data-testid="stDataFrame"]'); if (!df) return;
              const cv = df.querySelector('canvas'); if (!cv) return; const r = cv.getBoundingClientRect();
              const idx = Math.floor((e.clientY - r.top - 35) / 70); if (idx < 0) return;
              const inp = P.document.querySelector('.st-key-dbl_row input'); if (!inp) return;
              const set = Object.getOwnPropertyDescriptor(P.HTMLInputElement.prototype, 'value').set; inp.focus(); set.call(inp, idx + ':' + Date.now());
              inp.dispatchEvent(new P.Event('input', {bubbles: true})); inp.dispatchEvent(new P.KeyboardEvent('keypress', {key: 'Enter', code: 'Enter', keyCode: 13, charCode: 13, bubbles: true})); inp.blur(); }); }   // focus first: blur() on an unfocused input does nothing, so the value was never committed
        </script>""", height=0)

    # ====================================================================== OVERVIEW
    elif section == "Overview":
        lane, site, sig, topn, layout, _h = lens()
        site1 = site if site != "All sites" else pockets[0].pocket_id
        best = ranked(lane, site1, sig, topn)
        if best.empty: st.info("No poses."); st.stop()
        a, b = st.columns([1.15, 1]); a.altair_chart(top_bar(best, sig, lane), use_container_width=True); b.altair_chart(histogram(best, sig, lane), use_container_width=True)
        c, d = st.columns([1.15, 1]); c.altair_chart(size_scatter(best, sig, lane), use_container_width=True)
        wb = warhead_bar(res)
        if wb is not None: d.altair_chart(wb, use_container_width=True)
        else: d.info("No recognised warhead motifs in this library.")
        st.caption("The dashed line shows how the score tracks molecular weight; a steep slope means the ranking partly reflects ligand size (a known docking-score bias).")
        style.cards([("Ranked compounds", str(len(best)), site1), ("Top 10% by signal", str(int(best.top10.sum())), SIGNAL_LABEL[sig].split(' (')[0]),
                     ("Recognised warhead", str(int(best.warhead_classes.map(has_wh).sum())), "chemical annotation"),
                     ("Plausible covalent geometry", str(int(best.covalent_geometry.isin(['plausible', 'covalent_bond_modelled']).sum())), "pose geometry only")])

    # ====================================================================== POCKETS & CONFORMATIONS
    elif section == "Pockets & conformations":
        st.markdown("**All sites**")
        pr = pd.DataFrame([{"Site": p.pocket_id, "Source": p.source.replace("_", " "), "Status": "experimental ligand site" if p.experimentally_confirmed else "predicted / user-defined",
                            "Residues": len(p.residues), "Box (Å)": " × ".join(f"{x:.0f}" for x in p.size), "Compounds docked": int(docked[docked.pocket_id == p.pocket_id].ligand_key.nunique())} for p in pockets])
        st.dataframe(pr, hide_index=True, use_container_width=True)
        rec0 = R["receptor"]
        shared = [p for p in pockets if (p.receptor_path or rec0) == rec0]
        components.html(viewer3d.explorer_html([viewer3d.panel("Input structure", "all sites on the input conformation", rec0, [pocket_dict(p) for p in shared], [], target, None, None)], 600), height=760)
        lac = [p for p in pockets if str(p.metrics.get("tool", "")).startswith("Lacuna")]
        st.markdown("### Lacuna conformations")
        if not lac: st.info("Lacuna was not run in this screen. Select Lacuna under 'Predict pockets' on the New screen page to compare conformations.")
        else:
            lp = st.selectbox("Lacuna site", lac, format_func=lambda p: p.pocket_id)
            vb = {int(k): v for k, v in (lp.metrics.get("volume_by_conformer") or {}).items()}
            if vb:
                import altair as alt
                vd = pd.DataFrame({"conformer": list(vb), "volume": list(vb.values())}).sort_values("conformer")
                st.altair_chart(alt.Chart(vd).mark_bar(color=style.BLUE, size=14, cornerRadiusEnd=3).encode(x=alt.X("conformer:O", title="Conformer (0 = input structure)"), y=alt.Y("volume:Q", title="Pocket volume (Å³)"),
                                tooltip=["conformer", "volume"]).properties(height=260, title=f"{lp.pocket_id}: pocket volume across the ensemble (conformers where it appears)"), use_container_width=True)
            ens = lp.metrics.get("ensemble_pdb"); idxs = sorted(vb) or [lp.metrics.get("open_conformer_idx", 0)]
            ci = st.select_slider("Conformer to compare with the input structure", options=idxs, value=lp.metrics.get("open_conformer_idx", idxs[0])) if len(idxs) > 1 else idxs[0]
            if ens and Path(ens).exists():
                tmpf = rd / "renders"; tmpf.mkdir(exist_ok=True); cf = tmpf / f"conf_{lp.pocket_id}_{ci}.pdb"; cf.write_text(conformer_text(ens, int(ci)))
                A = viewer3d.panel("Input structure (closed)", ("pocket not detected in this conformation" if not lp.metrics.get('volume_input_A3') else f"pocket volume {lp.metrics.get('volume_input_A3')} Å³"), rec0, [pocket_dict(lp)], [], target, None, lp.pocket_id)
                B = viewer3d.panel(f"Lacuna conformation {ci}", f"pocket volume {vb.get(ci, '—')} Å³", str(cf), [pocket_dict(lp)], [], target, None, lp.pocket_id)
                components.html(viewer3d.explorer_html([A, B], 640), height=800)
            else: st.caption("Ensemble file not available for this run (re-run the screen to save per-conformer data).")
            style.note(lp.metrics.get("warning", ""), "note")

    # ====================================================================== EXPERIMENTAL
    elif section == "Experimental data":
        exp = R["experimental"]
        if exp is None: st.info("No experimental data was supplied. Add a CSV (compound_id + measured columns) on the New screen page.")
        else:
            st.dataframe(exp, hide_index=True, use_container_width=True)
            e1, e2, e3 = st.columns(3)
            lab = e1.selectbox("Label column", [c for c in exp.columns if c != "compound_id"])
            kind = e2.radio("Label type", ["binary (active = 1)", "continuous (higher = more potent)"]); ep = e3.selectbox("Site", [p.pocket_id for p in pockets])
            b = best_pose_per_compound(res, "gnina_cnn_vs", MODE_NONCOVALENT); b = b[b.pocket_id == ep]
            st.dataframe(evaluate(b, exp, lab, "binary" if kind.startswith("binary") else "continuous"), hide_index=True, use_container_width=True)
            style.note("Metrics describe agreement with the supplied labels on this set only. A signal that does not beat the <b>baseline_mol_weight</b> row is mostly measuring ligand size. Small sets give wide intervals.", "note")

    # ====================================================================== QC
    elif section == "Run status & QC":
        st.markdown("**Every input record, accounted for**")
        sc = res.drop_duplicates("ligand_key").groupby("processing_status").size().rename("compounds").reset_index()
        q1, q2 = st.columns([1, 2]); q1.dataframe(sc, hide_index=True, use_container_width=True)
        bad = res[res.processing_status != "docked"]
        q2.dataframe(bad[[c for c in ["compound_id", "processing_status", "failure_reason", "warnings"] if c in bad.columns]], hide_index=True, use_container_width=True)
        if "partial_failures" in res and res.partial_failures.notna().any():
            st.markdown("**Partial failures**"); st.dataframe(res[res.partial_failures.notna()][["compound_id", "partial_failures"]].drop_duplicates(), hide_index=True, use_container_width=True)
        if "covalent_note" in res and res.covalent_note.notna().any():
            st.markdown("**Covalent lane: no valid pose**"); st.dataframe(res[res.covalent_note.notna()][["compound_id", "covalent_note"]].drop_duplicates(), hide_index=True, use_container_width=True)
        if R["covalent_report"] is not None: st.markdown("**Covalent path notes**"); st.dataframe(R["covalent_report"], hide_index=True, use_container_width=True)
        if R["admet_meta"]:
            with st.expander("Property model metadata"): st.json(R["admet_meta"])
        with st.expander("Reproducibility manifest (versions, hashes, config)"): st.json(R["manifest"] or meta)

    # ====================================================================== EXPORT
    else:
        e1, e2 = st.columns(2)
        with e1:
            st.markdown("**Output package**")
            from pscreen import api
            try:
                od = rd / "outputs"
                if not (od / "REPORT.md").exists(): api.export_outputs(rd)
                zb = io.BytesIO()
                with zipfile.ZipFile(zb, "w", zipfile.ZIP_DEFLATED) as zf:
                    for f in sorted(od.iterdir()): zf.write(f, f"{rd.name}/{f.name}")
                st.download_button("⬇ Everything (REPORT.md, ranked CSV, status, sites, best poses SDF)", zb.getvalue(), f"{rd.name}_outputs.zip", mime="application/zip", type="primary")
            except Exception as ex: st.caption(f"Output package unavailable: {ex}")
            st.markdown("**Tables & records**")
            for f, lbl in [("results_all_poses.csv", "All poses (CSV)"), ("results_best_pose.csv", "Best pose per compound / site / lane (CSV)"), ("run_config.json", "Reproducible run config (JSON)"), ("run_manifest.json", "Manifest with tool versions (JSON)")]:
                if (rd / f).exists(): st.download_button(lbl, (rd / f).read_bytes(), f)
        with e2:
            st.markdown("**Structures**")
            z = rd / "structures_bundle.zip"
            if z.exists(): st.download_button("Bundle: all poses + receptor(s) + sites + config (zip)", z.read_bytes(), f"{rd.name}_structures.zip")
            sel = st.multiselect("Pick compounds for one SDF of their best poses", sorted(docked.compound_id.unique()))
            if sel:
                buf = io.StringIO(); w = Chem.SDWriter(buf)
                for _, r in best_pose_per_compound(res, "gnina_cnn_vs", None).query("compound_id in @sel").iterrows():
                    m = pose_block(r.pose_sdf, int(r.pose_idx)); m.SetProp("_Name", f"{r.compound_id}|{r.pocket_id}|{r.docking_mode}")
                    for k in ["gnina_minimized_affinity_kcal", "gnina_cnn_pose_score", "gnina_cnn_affinity_pK", "docking_mode", "pocket_id", "pocket_evidence"]: m.SetProp(k, str(r[k]))
                    w.write(m)
                w.close(); st.download_button("Download selected poses (SDF)", buf.getvalue(), "selected_poses.sdf")
            st.caption("PyMOL session (.pse) / script (.pml) downloads are under 'Other viewers' on the Explore page.")
