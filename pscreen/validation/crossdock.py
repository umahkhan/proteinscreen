"""Cross-target enrichment: each target's own crystal ligand versus the ligands of the other targets (all docked into every pocket).

STAND-IN design: DUD-E (blocked to scripted download), PDBbind (registration) and Binding MOAD (retired) were not available. It is weaker
than a property-matched decoy benchmark: (a) the 'decoys' are not guaranteed non-binders, (b) they are not property-matched, so docking-score
size bias is NOT controlled for -> a heavy-atom-count baseline is reported next to every signal, (c) there is one positive per target.
Metric per target: rank of the native ligand among all ligands (1 = best); summary = native-ligand top-1/top-3 rate and mean rank percentile.
Chance top-1 = 1/N. Scores are only ever compared within one pocket, never across pockets.
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
from ..docking import DockConfig, run_screen, Store
from ..runio import load_pockets

SIGNALS = {"gnina_minimized_affinity_kcal": ("minimized_affinity", True), "gnina_cnn_pose_score": ("cnn_score", False),
           "gnina_cnn_affinity_pK": ("cnn_affinity", False), "gnina_cnn_vs": ("cnn_vs", False)}


def run(redock_dir: str | Path, cfg: DockConfig, progress=None) -> pd.DataFrame:
    rd = Path(redock_dir).resolve()
    store = Store(rd / "results.sqlite")
    meta = {k: json.loads(v) for k, v in store.db.execute("SELECT key,value FROM meta")}
    pockets = load_pockets(meta["pockets"]); prepared = pd.read_csv(rd / "prepared.csv")
    run_screen(prepared, meta["receptor"], pockets, cfg, rd, progress=progress)         # resumable: native-pair jobs already done are skipped
    return score(rd, prepared, pockets)


def score(rd: Path, prepared: pd.DataFrame, pockets) -> pd.DataFrame:
    store = Store(rd / "results.sqlite"); poses, jobs = store.poses(), store.jobs()
    heavy = prepared.set_index("ligand_key").heavy_atoms
    rows = []
    for pk in pockets:
        native = pk.pocket_id.replace("_site", "")
        g = poses[(poses.pocket_id == pk.pocket_id) & (poses["mode"] == "noncovalent")]
        n_failed = int(((jobs.pocket_id == pk.pocket_id) & (jobs.status == "failed")).sum())
        for name, (col, asc) in SIGNALS.items():
            best = g.sort_values(col, ascending=asc).groupby("ligand_key").head(1).set_index("ligand_key")[col]
            if native not in best.index: rows.append({"pocket": native, "signal": name, "note": "native ligand has no pose"}); continue
            order = best.sort_values(ascending=asc).index.tolist()
            rows.append({"pocket": native, "signal": name, "n_ligands": len(best), "native_rank": order.index(native) + 1, "n_failed_jobs": n_failed})
        hb = heavy.reindex(g.ligand_key.unique()).sort_values(ascending=False)       # baseline: bigger ligand ranks higher
        if native in hb.index:
            rows.append({"pocket": native, "signal": "baseline_heavy_atoms", "n_ligands": len(hb), "native_rank": int(hb.index.tolist().index(native) + 1), "n_failed_jobs": n_failed})
    res = pd.DataFrame(rows); res.to_csv(rd / "crossdock_ranks.csv", index=False)
    return res


def summarize(res: pd.DataFrame) -> pd.DataFrame:
    r = res.dropna(subset=["native_rank"])
    out = r.groupby("signal").apply(lambda d: pd.Series({"n_targets": len(d), "top1": int((d.native_rank == 1).sum()), "top3": int((d.native_rank <= 3).sum()),
                                                         "mean_rank_percentile": round(float(((d.native_rank - 1) / (d.n_ligands - 1)).mean()), 3),
                                                         "chance_top1": round(float((1 / d.n_ligands).mean()), 3)}), include_groups=False)
    return out
