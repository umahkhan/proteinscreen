"""Product-facing layer: validate inputs, create a run, and write the standard output package.

INPUTS   protein (PDB/mmCIF) + compounds (.csv/.tsv/.smi/.sdf) [+ target residue(s)] [+ known ligand | box] [+ pocket predictors] [+ ADMET]
OUTPUTS  runs/<name>/outputs/
           REPORT.md              human-readable summary: inputs, sites, top compounds, failures, limitations
           summary.json           the same facts, machine-readable
           ranked_compounds.csv   one row per compound and docking lane, ranked by GNINA CNN VS (never mixes covalent / non-covalent)
           compound_status.csv    EVERY input compound with its status (docked / invalid / duplicate / failed) and reason
           pockets.csv            binding sites used (tool, residues, centre) and the pocket-tool status
           best_poses.sdf         best pose per compound and lane (3D, properties attached)
The UI and the CLI both go through here, so they behave identically.
"""
from __future__ import annotations
import json, re, tempfile
from pathlib import Path
import pandas as pd
from rdkit import Chem

from .pipeline import RunConfig, prepare_receptor
from .prep import read_library
from .residues import validate_target
from .results import best_pose_per_compound

LIB_EXT = (".csv", ".tsv", ".smi", ".txt", ".smiles", ".sdf", ".sd")
MAX_LIBRARY = 1000
LIMITS = ("Docking scores are computational predictions; they do not establish binding, covalent engagement, potency or selectivity. "
          "Covalent and non-covalent scores are on different scales and are not compared. Pockets from apo structures and cryptic-pocket predictors are hypotheses.")


def _s(v) -> str:
    return "" if v is None or (isinstance(v, float) and v != v) else str(v)


def slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")[:60]


def preflight(protein: str, compounds: str, residues: list[str] | None = None, known_ligand: str = "", box_center: list | None = None) -> dict:
    """Cheap checks before any compute. Returns {ok, errors[], warnings[], info{}}; never raises for bad user input."""
    errors, warnings, info = [], [], {}
    p, c = Path(protein), Path(compounds)
    if not p.exists(): errors.append(f"Protein file not found: {protein}")
    if not c.exists(): errors.append(f"Compound file not found: {compounds}")
    elif c.suffix.lower() not in LIB_EXT: errors.append(f"Unsupported compound format '{c.suffix}' (use {', '.join(LIB_EXT)})")
    if known_ligand and not Path(known_ligand).exists(): errors.append(f"Known-ligand file not found: {known_ligand}")
    if box_center and len(box_center) != 3: errors.append("box_center needs three numbers x,y,z")
    rec = None
    if p.exists():
        try:
            tmp = Path(tempfile.mkdtemp()) / "receptor.pdb"; info["protein"] = prepare_receptor(p, tmp); rec = tmp
            if info["protein"]["n_residues"] < 30: warnings.append(f"Only {info['protein']['n_residues']} residues: is this a full protein?")
        except Exception as e: errors.append(str(e))
    if c.exists() and c.suffix.lower() in LIB_EXT:
        try:
            recs = read_library(c); bad = [r for r in recs if r.error]
            info["compounds"] = {"records": len(recs), "unparseable": len(bad)}
            if not recs: errors.append("The compound file has no records")
            elif len(bad) == len(recs): errors.append("None of the compounds could be parsed: " + bad[0].error)
            elif bad: warnings.append(f"{len(bad)} of {len(recs)} compounds cannot be parsed and will be reported as invalid (e.g. {bad[0].original_id}: {bad[0].error})")
            if len(recs) > MAX_LIBRARY: warnings.append(f"{len(recs)} compounds (> {MAX_LIBRARY}): the run will be slow")
        except Exception as e: errors.append(str(e))
    if rec is not None:
        for r in residues or []:
            v = validate_target(str(rec), r)
            if not v["ok"]: errors.append(f"Target residue '{r}': " + " ".join(v["errors"]))
            else: warnings += v["warnings"]
    return {"ok": not errors, "errors": errors, "warnings": warnings, "info": info}


def create_run(name: str, protein: str, compounds: str, root: str | Path = "runs", residues: list[str] | None = None, known_ligand: str = "",
               box_center: list | None = None, box_size: list | None = None, pockets: list[str] | None = None, admet: bool = False,
               covalent: bool = True, fast: bool = False, **extra) -> Path:
    """Write runs/<slug>/run_config.json (input paths made absolute). Run it with pscreen.pipeline.run(rundir) or the UI launcher."""
    key = slug(name)
    if not key: raise ValueError("Run name is empty")
    rd = Path(root) / key
    if (rd / "run_config.json").exists(): raise FileExistsError(f"Run '{key}' already exists: pick another name or delete {rd}")
    rd.mkdir(parents=True, exist_ok=True)
    res = list(residues or [])
    cfg = RunConfig(name=name, receptor=str(Path(protein).resolve()), library=str(Path(compounds).resolve()), target_residue=res[0] if res else "", target_residues=res[1:],
                    known_ligand=str(Path(known_ligand).resolve()) if known_ligand else "", box_center=list(box_center or []), covalent=covalent and bool(res),
                    predict_pockets=bool(pockets), predictors=list(pockets or ["p2rank"]), chemprop=admet)
    if box_size: cfg.box_size = list(box_size)
    if fast: cfg.exhaustiveness, cfg.num_modes = 4, 3
    for k, v in extra.items():
        if k not in RunConfig.__dataclass_fields__: raise TypeError(f"unknown option {k}")
        setattr(cfg, k, v)
    cfg.save(rd); return rd


# ------------------------------------------------------------------------------------------------------------------ outputs
def ranked_compounds(res: pd.DataFrame, signal: str = "gnina_cnn_vs") -> pd.DataFrame:
    """Best pose per compound within each docking lane (across all sites), ranked inside the lane."""
    b = best_pose_per_compound(res, signal, None)
    if b.empty: return b
    b = b.sort_values(signal, ascending=signal == "gnina_minimized_affinity_kcal").drop_duplicates(["docking_mode", "ligand_key"])
    b["rank"] = b.groupby("docking_mode")[signal].rank(ascending=signal == "gnina_minimized_affinity_kcal", method="first").astype(int)
    cols = ["docking_mode", "rank", "compound_id", "pocket_id", "pocket_source", "gnina_cnn_vs", "gnina_cnn_pose_score", "gnina_cnn_affinity_pK", "gnina_minimized_affinity_kcal",
            "warhead_classes", "target_residue", "reactive_to_target_distance", "covalent_geometry", "mol_weight", "clogp", "tpsa", "canonical_smiles", "pose_sdf"]
    cols += [c for c in b.columns if c.startswith("chemprop__")] + [c for c in b.columns if c.startswith("exp_")]
    return b[[c for c in cols if c in b.columns]].sort_values(["docking_mode", "rank"]).reset_index(drop=True)


def export_outputs(rundir: str | Path) -> Path:
    """Write the standard output package (see module docstring). Safe to call again; works on finished or partial runs."""
    from .runio import load_run
    rd = Path(rundir); R = load_run(rd); res = R["results"]; out = rd / "outputs"; out.mkdir(exist_ok=True)
    rk = ranked_compounds(res); rk.drop(columns=["pose_sdf"], errors="ignore").to_csv(out / "ranked_compounds.csv", index=False)
    first = res.drop_duplicates("ligand_key")
    stat = (R["prepared"][["compound_id", "status", "reason", "canonical_smiles", "warhead_classes"]]
            .rename(columns={"status": "prep_status", "reason": "prep_reason"}))
    dock = res.groupby("compound_id").processing_status.agg(lambda s: "docked" if (s == "docked").any() else s.iloc[0]).rename("final_status")
    why = first.set_index("compound_id")[["failure_reason"]] if "failure_reason" in first else None
    stat = stat.merge(dock, on="compound_id", how="left")
    if why is not None: stat = stat.merge(why, on="compound_id", how="left")
    stat.to_csv(out / "compound_status.csv", index=False)
    pk = pd.DataFrame([{"pocket_id": p.pocket_id, "source": p.source, "tool": p.metrics.get("tool", ""), "n_residues": len(p.residues), "center_x": p.center[0], "center_y": p.center[1],
                        "center_z": p.center[2], "residues": " ".join(p.residues)} for p in R["pockets"]])
    pk.to_csv(out / "pockets.csv", index=False)
    n = 0
    with Chem.SDWriter(str(out / "best_poses.sdf")) as w:
        for _, r in rk.iterrows():
            ms = [m for m in Chem.SDMolSupplier(r.pose_sdf, removeHs=False)] if Path(str(r.pose_sdf)).exists() else []
            full = res[(res.compound_id == r.compound_id) & (res.docking_mode == r.docking_mode) & (res.pocket_id == r.pocket_id)].sort_values("gnina_cnn_vs", ascending=False)
            if not ms or full.empty: continue
            m = ms[int(full.iloc[0].pose_idx)]; m.SetProp("_Name", str(r.compound_id))
            for k in ("docking_mode", "rank", "pocket_id", "gnina_cnn_vs", "gnina_minimized_affinity_kcal", "warhead_classes", "covalent_geometry"):
                if k in r and pd.notna(r[k]): m.SetProp(k, str(r[k]))
            w.write(m); n += 1
    jobs = R["store"].jobs(); mf = R["manifest"]; cfg = R["config"]
    summary = {"name": cfg.get("name"), "finished": mf.get("finished"), "protein": Path(cfg["receptor"]).name, "protein_info": mf.get("receptor_info"), "compounds_file": Path(cfg["library"]).name,
               "library": mf.get("library_summary"), "target_residues": [r for r in [cfg.get("target_residue"), *cfg.get("target_residues", [])] if r],
               "pocket_tools": mf.get("pocket_tools"), "n_pockets": len(R["pockets"]), "jobs": {f"{m}/{s}": int(c) for (m, s), c in jobs.groupby(["mode", "status"]).size().items()},
               "n_ranked": {m: int((rk.docking_mode == m).sum()) for m in rk.docking_mode.unique()} if len(rk) else {}, "poses_written": n,
               "failed_jobs": jobs[jobs.status == "failed"][["ligand_key", "pocket_id", "mode", "error"]].to_dict("records"), "gnina": mf.get("gnina"), "covalent_error": mf.get("covalent_error"),
               "limitations": LIMITS}
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    L = [f"# pscreen report — {summary['name']}", "", f"Finished {summary['finished']}  ·  protein `{summary['protein']}`  ·  compounds `{summary['compounds_file']}`  ·  {summary['gnina']}", "",
         "## Inputs", f"- Library: {(summary['library'] or {}).get('total', '?')} compounds, {(summary['library'] or {}).get('ok', '?')} usable, {(summary['library'] or {}).get('with_warhead', 0)} with a recognised warhead", f"- Target residue(s): {', '.join(summary['target_residues']) or '—'}", "", "## Binding sites"]
    L += [f"- `{k}`: {v}" for k, v in (summary["pocket_tools"] or {}).items()] + [f"- {len(R['pockets'])} site(s) docked (see pockets.csv)", "", "## Top compounds (per docking lane)"]
    for m in rk.docking_mode.unique() if len(rk) else []:
        t = rk[rk.docking_mode == m].head(10)
        L += ["", f"### {m}", "", "| rank | compound | site | CNN VS | affinity kcal | warhead | covalent geometry |", "|---|---|---|---|---|---|---|"]
        L += [f"| {r['rank']} | {r.compound_id} | {r.pocket_id} | {r.gnina_cnn_vs:.2f} | {r.gnina_minimized_affinity_kcal:.2f} | {_s(r.get('warhead_classes'))} | {_s(r.get('covalent_geometry'))} |" for _, r in t.iterrows()]
    if not len(rk): L += ["", "_No compound produced a docked pose._"]
    bad = stat[stat.final_status != "docked"]
    L += ["", "## Not docked / problems", ""] + ([f"- {r.compound_id}: {r.final_status} — {_s(r.get('failure_reason')) or _s(r.prep_reason)}" for _, r in bad.iterrows()] or ["- none"])
    if summary["failed_jobs"]: L += [f"- failed docking jobs: {len(summary['failed_jobs'])} (see summary.json)"]
    if summary["covalent_error"]: L += [f"- covalent stage error: {summary['covalent_error']}"]
    L += ["", "## Limitations", "", LIMITS, ""]
    (out / "REPORT.md").write_text("\n".join(L))
    return out
