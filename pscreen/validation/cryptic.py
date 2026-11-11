"""Cryptic-pocket discovery evaluation on the PocketMiner dataset (Meller et al. 2023; bundled in bowman-lab/ae-pocketminer, MIT).

Input is the cleaned APO structure only. Ground truth = residues labelled 1 (they form the pocket in a corresponding holo structure);
2 = ambiguous (excluded); 0 = confident negative (negative proteins). Labels are indexed by protein residue order in the structure.
Pocket-level hit (Lacuna's size-robust criterion): residue Jaccard with the true site >= 0.25 OR centroid within 4 A of the true site's
heavy-atom centroid. Reported for top-1/3/5 pockets per tool. Pocket DISCOVERY is evaluated here, separately from compound ranking.
Caveats: PocketMiner's own test split is held out for PocketMiner; Lacuna's ranker was fitted on CryptoBench folds and P2Rank on its own
sets; overlap of these benchmarks with those training sets is not verified. The PocketMiner threshold (default 0.7) was an a-priori choice;
scripts_cryptic.py selects it on the VAL split only (never on test). Small n (24 positive proteins across val+test): wide uncertainty.
"""
from __future__ import annotations
import glob
import os
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from ..pocket_tools import run_p2rank, run_pocketminer, run_lacuna, _residue_table
from ..docking import Pocket

DS = Path(__file__).resolve().parents[2] / "third_party" / "ae-pocketminer" / "data" / "pm-dataset"


def load_split(split: str) -> list[dict]:
    ids = np.load(DS / f"{split}_apo_ids_with_chainids.npy", allow_pickle=True)
    labs = np.load(DS / f"{split}_label_dictionary.npy", allow_pickle=True).item()
    files = {os.path.basename(f).lower(): f for f in glob.glob(str(DS / "apo-structures" / "*.pdb"))}
    out = []
    for i in map(str, ids):
        f = files.get(f"{i.lower()}_clean_h.pdb")
        key = next((k for k in (i, i[:4], i.lower(), i[:4].lower()) if k in labs), None)
        out.append({"id": i, "file": f, "labels": None if key is None else np.array(labs[key]).astype(int),
                    "error": None if (f and key) else ("apo structure file not found" if not f else "label key not found in label dictionary")})
    return out


def _sites(entry: dict, res_rows: pd.DataFrame, table: dict) -> tuple[set, np.ndarray | None]:
    lab = entry["labels"]
    pos = {(r.chain, int(r.resseq)) for r, l in zip(res_rows.itertuples(), lab) if l == 1}
    pos = {k for k in pos if k in table}
    if not pos: return pos, None
    return pos, np.vstack([table[k]["xyz"] for k in pos]).mean(0)


def _parse(res: str) -> tuple[str, int]:
    ch, rest = res.split(":"); return ch, int("".join(c for c in rest if c.isdigit() or c == "-"))


def hit(pocket: Pocket, true: set, tcenter: np.ndarray) -> tuple[bool, float, float]:
    ps = {_parse(r) for r in pocket.residues}
    jac = len(ps & true) / max(1, len(ps | true)); d = float(np.linalg.norm(np.array(pocket.center) - tcenter))
    return (jac >= 0.25 or d <= 4.0), round(jac, 3), round(d, 2)


def run(split: str, outdir: str | Path, threshold: float = 0.7, tools=("p2rank", "pocketminer", "lacuna"), progress=None) -> dict:
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    pocket_rows, resid_rows, status = [], [], []
    for e in load_split(split):
        pid = e["id"]
        if e["error"]: status.append({"id": pid, "status": "skipped", "reason": e["error"]}); continue
        work = outdir / pid; table = _residue_table(e["file"]); lab = e["labels"]
        is_pos = bool((lab == 1).any())
        row_status = {"id": pid, "status": "ok", "reason": None, "kind": "positive" if is_pos else "negative"}
        try:
            pm_pockets, csvp = run_pocketminer(e["file"], work / "pm", threshold=threshold, max_pockets=5)
            res = pd.read_csv(csvp)
            if len(res) != len(lab): raise ValueError(f"label length {len(lab)} != residues {len(res)}")
            true, tc = _sites(e, res, table)
            for r, l in zip(res.itertuples(), lab):
                if l in (0, 1): resid_rows.append({"id": pid, "label": int(l), "probability": r.probability})
            row_status.update(n_true_residues=len(true), max_prob=float(res.probability.max()), n_pm_pockets=len(pm_pockets))
            if is_pos and tc is not None:
                allp = {"pocketminer": pm_pockets}
                if "p2rank" in tools:
                    try: allp["p2rank"] = run_p2rank(e["file"], work / "p2rank", max_pockets=5)
                    except Exception as ex: row_status["p2rank_error"] = str(ex)[:120]
                if "lacuna" in tools:
                    try: allp["lacuna"] = run_lacuna(e["file"], work / "lacuna", n_conformers=10, top=5)
                    except Exception as ex: row_status["lacuna_error"] = str(ex)[:120]
                for tool, pks in allp.items():
                    for rank, pk in enumerate(pks, 1):
                        h, jac, d = hit(pk, true, tc)
                        pocket_rows.append({"id": pid, "tool": tool, "rank": rank, "hit": h, "jaccard": jac, "centroid_dist_A": d, "n_residues": len(pk.residues)})
        except Exception as ex:
            row_status.update(status="failed", reason=str(ex)[:200])
        status.append(row_status)
        if progress: progress(pid, row_status["status"])
    pd.DataFrame(pocket_rows).to_csv(outdir / "pockets.csv", index=False)
    pd.DataFrame(resid_rows).to_csv(outdir / "residues.csv", index=False)
    pd.DataFrame(status).to_csv(outdir / "status.csv", index=False)
    return {"pockets": pd.DataFrame(pocket_rows), "residues": pd.DataFrame(resid_rows), "status": pd.DataFrame(status)}


def summarize(r: dict, thr_neg: float = 0.7) -> dict:
    st, pk, rs = r["status"], r["pockets"], r["residues"]
    ok = st[st.status == "ok"]
    out = {"n_proteins": len(st), "n_ok": len(ok), "n_skipped_or_failed": int((st.status != "ok").sum()),
           "n_positive_proteins": int((ok.kind == "positive").sum()), "n_negative_proteins": int((ok.kind == "negative").sum())}
    if len(rs) and rs.label.nunique() == 2:
        out["pocketminer_residue_roc_auc_pooled"] = round(float(roc_auc_score(rs.label, rs.probability)), 3)
        out["pocketminer_residue_pr_auc_pooled"] = round(float(average_precision_score(rs.label, rs.probability)), 3)
    neg = ok[ok.kind == "negative"]
    if len(neg): out["pocketminer_negative_proteins_with_any_residue_above_threshold"] = f"{int((neg.max_prob >= thr_neg).sum())}/{len(neg)}"
    if len(pk):
        for tool, g in pk.groupby("tool"):
            n = g.id.nunique()
            for k in (1, 3, 5):
                out[f"{tool}_top{k}_recovery"] = round(float(g[g["rank"] <= k].groupby("id").hit.any().reindex(g.id.unique(), fill_value=False).mean()), 3)
            out[f"{tool}_n_proteins_with_pockets"] = n
    return out

