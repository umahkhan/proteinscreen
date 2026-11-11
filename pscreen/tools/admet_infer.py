"""Runs in the isolated ADMET-AI env (admet-ai 2.0.1 / Chemprop 2.x, MIT). Usage: admet_infer.py <in.csv(ligand_key,smiles)> <out.csv> <endpoint,...>
Writes one column per endpoint (classification -> probability, regression -> value in the endpoint's units)."""
import sys, warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
import pandas as pd
from admet_ai import ADMETModel
inp, out, eps = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
df = pd.read_csv(inp)
preds = ADMETModel(num_workers=0).predict(smiles=df["smiles"].tolist())
missing = [e for e in eps if e not in preds.columns]
if missing: sys.exit(f"endpoints not available: {missing}")
res = preds[eps].reset_index(drop=True); res.insert(0, "ligand_key", df["ligand_key"].values)
res.to_csv(out, index=False)
