"""Retrospective evaluation of ranking signals against user-supplied experimental data.

Signals are evaluated one at a time (no composite). A molecular-weight / heavy-atom baseline is always included because docking
scores correlate with ligand size; a signal that does not beat that baseline is mostly measuring size. Orientation: for
'gnina_minimized_affinity_kcal' lower is better (negated internally); for all other signals higher is better.
Binary labels -> ROC-AUC, PR-AUC, enrichment factors, top-k precision.  Continuous labels (e.g. pIC50, higher = more potent) -> Spearman.
These metrics describe agreement with the supplied labels on THIS set; with few compounds the uncertainty is large (bootstrap CI reported).
"""
from __future__ import annotations
import warnings
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score

SIGNALS = {"gnina_minimized_affinity_kcal": -1, "gnina_cnn_pose_score": 1, "gnina_cnn_affinity_pK": 1, "gnina_cnn_vs": 1,
           "baseline_mol_weight": 1, "baseline_heavy_atoms": 1}


def _ef(y: np.ndarray, s: np.ndarray, frac: float) -> float | None:
    n = len(y); k = max(1, int(round(n * frac)))
    if y.sum() == 0: return None
    top = np.argsort(-s, kind="stable")[:k]
    return float((y[top].mean()) / (y.mean()))


def _boot(fn, y, s, n=500, seed=0):
    rng = np.random.default_rng(seed); vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(set(y[i])) < 2: continue            # resample with one class: AUC undefined, skip silently
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try: vals.append(fn(y[i], s[i]))
            except Exception: pass
    return (round(float(np.percentile(vals, 2.5)), 3), round(float(np.percentile(vals, 97.5)), 3)) if len(vals) > 20 else (None, None)


def evaluate(best: pd.DataFrame, exp: pd.DataFrame, label_col: str, kind: str = "binary", signals: list[str] | None = None) -> pd.DataFrame:
    """best: one row per compound (e.g. best_pose_per_compound for ONE pocket and ONE docking mode); exp: compound_id + label_col."""
    d = best.merge(exp[["compound_id", label_col]], on="compound_id", how="inner").dropna(subset=[label_col])
    for src, dst in (("mol_weight", "baseline_mol_weight"), ("heavy_atoms", "baseline_heavy_atoms")):
        if src in d.columns and d[src].notna().any(): d[dst] = d[src]
    rows = []
    y = d[label_col].astype(float).values
    for sig in signals or [s for s in SIGNALS if s in d.columns]:
        s = SIGNALS.get(sig, 1) * d[sig].astype(float).values
        ok = ~np.isnan(s)
        r = {"signal": sig, "n": int(ok.sum())}
        if ok.sum() < 3: r.update(note="fewer than 3 compounds with this signal"); rows.append(r); continue
        if kind == "binary":
            yy = y[ok].astype(int); ss = s[ok]
            if yy.min() == yy.max(): r.update(note="only one class present"); rows.append(r); continue
            r.update(roc_auc=round(float(roc_auc_score(yy, ss)), 3), pr_auc=round(float(average_precision_score(yy, ss)), 3),
                     base_rate=round(float(yy.mean()), 3), ef1pct=_ef(yy, ss, 0.01), ef5pct=_ef(yy, ss, 0.05))
            r["roc_auc_ci95"] = _boot(roc_auc_score, yy, ss)
        else:
            rho = spearmanr(y[ok], s[ok]); r.update(spearman=round(float(rho.statistic), 3), p=float(rho.pvalue))
        rows.append(r)
    return pd.DataFrame(rows)
