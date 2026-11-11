"""Optional ADMET property predictions via ADMET-AI (Chemprop-RDKit ensembles pretrained on Therapeutics Data Commons datasets).

What this is: pretrained, third-party ADMET models (no training here). What this is NOT: a predictor of binding affinity, covalent
reactivity, potency or selectivity. No generic Chemprop checkpoint provides those, and none is used. These columns are kept separate
from GNINA scores and never combined with them.
Endpoints with negative held-out R^2 in ADMET-AI's own evaluation (e.g. half-life, VDss) are excluded by default.
If the isolated env is missing, predictions are disabled with a clear message and screening continues.
"""
from __future__ import annotations
import json, subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PY = ROOT / "envs" / "admet" / "bin" / "python"
INFER = ROOT / "pscreen" / "tools" / "admet_infer.py"
DEFAULT_ENDPOINTS = ["Solubility_AqSolDB", "Lipophilicity_AstraZeneca", "Caco2_Wang", "PPBR_AZ", "hERG", "AMES", "DILI", "CYP3A4_Veith"]
MODEL_NAME = "ADMET-AI (Chemprop-RDKit D-MPNN ensemble, 5 models/endpoint)"
LIMITATIONS = ("Trained on public TDC datasets (see ADMET-AI paper, Bioinformatics 2024); applicability domain is drug-like small molecules. "
               "Does not predict target binding, covalent reactivity, potency or selectivity. Classification outputs are probabilities, "
               "not calibrated risk. Overlap between TDC training data and your compounds is not checked.")


def available() -> bool:
    return PY.exists()


def _version() -> str:
    return subprocess.run([str(PY), "-c", "import importlib.metadata as m;print(m.version('admet-ai'), 'chemprop', m.version('chemprop'))"],
                          capture_output=True, text=True).stdout.strip()


def _endpoint_table() -> pd.DataFrame:
    p = subprocess.run([str(PY), "-c", "import admet_ai,os;print(os.path.join(os.path.dirname(admet_ai.__file__),'resources','data','admet.csv'))"],
                       capture_output=True, text=True).stdout.strip()
    return pd.read_csv(p).set_index("id")


def predict(prepared: pd.DataFrame, outdir: str | Path, endpoints: list[str] | None = None) -> tuple[pd.DataFrame | None, dict]:
    """Returns (predictions keyed by ligand_key with columns chemprop__<endpoint>, metadata dict). (None, {disabled: reason}) if unavailable."""
    if not available():
        return None, {"disabled": f"ADMET-AI environment not found ({PY}); property prediction disabled. See docs/SETUP.md."}
    endpoints = endpoints or DEFAULT_ENDPOINTS
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    ok = prepared[prepared.status == "ok"][["ligand_key", "canonical_smiles"]].rename(columns={"canonical_smiles": "smiles"})
    ok.to_csv(outdir / "admet_in.csv", index=False)
    p = subprocess.run([str(PY), str(INFER), str(outdir / "admet_in.csv"), str(outdir / "admet_out.csv"), ",".join(endpoints)],
                       capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        return None, {"disabled": "ADMET-AI failed: " + "\n".join((p.stderr or p.stdout).strip().splitlines()[-3:])}
    res = pd.read_csv(outdir / "admet_out.csv").rename(columns={e: f"chemprop__{e}" for e in endpoints})
    tab = _endpoint_table()
    meta = {"model": MODEL_NAME, "version": _version(), "limitations": LIMITATIONS, "license": "MIT (admet-ai, chemprop)",
            "endpoints": {e: {"name": tab.loc[e, "name"], "task": tab.loc[e, "task_type"], "units": tab.loc[e, "units"] if tab.loc[e, "task_type"] == "regression" else "probability (0-1)",
                              "train_set_size_TDC": int(tab.loc[e, "size"]),
                              "heldout_test_metric": ({"AUROC": round(float(tab.loc[e, "AUROC"]), 3)} if tab.loc[e, "task_type"] == "classification" else {"R2": round(float(tab.loc[e, "R^2"]), 3), "MAE": round(float(tab.loc[e, "MAE"]), 3)})}
                          for e in endpoints}}
    (outdir / "admet_models.json").write_text(json.dumps(meta, indent=1))
    return res, meta
