"""Compact top-N compound table (images + key ADMET) and the full ADMET detail used in the compound popup."""
from __future__ import annotations
import json
import pandas as pd
from .depict import compound_png, warhead_png

ADMET_COMPACT = [("Solubility_AqSolDB", "Solubility logS", "num"), ("hERG", "hERG risk", "prob"), ("AMES", "AMES", "prob"), ("DILI", "DILI", "prob")]
BANDS = [(0.3, "low"), (0.7, "moderate"), (1.01, "high")]


def has_wh(v) -> bool: return isinstance(v, str) and v.strip() not in ("", "nan", "None")


def hits_of(prow) -> list[dict]:
    return json.loads(prow.warheads) if isinstance(prow.get("warheads"), str) else []


def build_table(R: dict, best: pd.DataFrame, top_n: int, show_site: bool) -> pd.DataFrame:
    prep = R["prepared"].set_index("ligand_key", drop=False)
    t = best.head(top_n).copy(); rows = []
    for _, r in t.iterrows():
        p = prep.loc[r.ligand_key]; hits = hits_of(p)
        m = tuple(sorted({a for h in hits for a in h["match_atoms"]})); rx = tuple(sorted({a for h in hits for a in h["reactive_atoms"]}))
        row = {"Rank": int(r["rank"]), "Compound": compound_png(p.ligand_sdf, m, rx, 130, 84), "Name": r.compound_id,
               "Warhead structure": warhead_png(p.ligand_sdf, tuple(hits[0]["match_atoms"]), tuple(hits[0]["reactive_atoms"]), 100, 70) if hits else None,
               "Warhead": hits[0]["warhead"].replace("_", " ") if hits else "—", "Site": r.pocket_id,
               "CNN pose": r.gnina_cnn_pose_score, "Vina-like": r.gnina_minimized_affinity_kcal, "CNN pK": r.gnina_cnn_affinity_pK,
               "Reactive→target (Å)": r.reactive_to_target_distance, "_pocket": r.pocket_id, "_ligand_key": r.ligand_key}
        for ep, label, _ in ADMET_COMPACT:
            c = f"chemprop__{ep}"
            if c in r.index: row[label] = r[c]
        rows.append(row)
    df = pd.DataFrame(rows)
    if not show_site and "Site" in df: df = df.drop(columns=["Site"])
    return df


def admet_rows(row, meta: dict | None) -> pd.DataFrame:
    eps = (meta or {}).get("endpoints", {})
    out = []
    for c in [c for c in row.index if str(c).startswith("chemprop__") and pd.notna(row[c])]:
        e = c.replace("chemprop__", ""); m = eps.get(e, {}); v = float(row[c]); cls = m.get("task") == "classification"
        out.append({"Endpoint": m.get("name", e), "Prediction": v, "Reads as": ("probability · " + next(b for lim, b in BANDS if v < lim)) if cls else m.get("units", ""),
                    "Model held-out performance": ", ".join(f"{k} {x}" for k, x in m.get("heldout_test_metric", {}).items()), "Training set (TDC)": m.get("train_set_size_TDC")})
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ RDKit chemistry columns and grid
EXTRA_COLUMNS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "Rot. bonds", "Lipinski violations", "Scaffold", "Similarity to selected"]


def add_chem_columns(R: dict, tab: pd.DataFrame, extra: list[str], primary_ligand_key: str | None) -> pd.DataFrame:
    """Add RDKit-derived columns on the prepared ligands: descriptors (Descriptors / Lipinski), rule-of-five violations, Bemis-Murcko scaffold,
    and Morgan(2) Tanimoto similarity to the selected compound."""
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    if not extra: return tab
    prep = R["prepared"].set_index("ligand_key", drop=False)
    mols = {k: Chem.RemoveHs(Chem.SDMolSupplier(prep.loc[k].ligand_sdf, removeHs=False)[0]) for k in tab["_ligand_key"].unique()}
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    fps = {k: fpg.GetFingerprint(m) for k, m in mols.items()} if "Similarity to selected" in extra else {}
    ref = fps.get(primary_ligand_key) if primary_ligand_key else None
    out = tab.copy()
    def col(name, fn): out[name] = [fn(mols[k]) for k in out["_ligand_key"]]
    if "MW" in extra: col("MW", lambda m: round(Descriptors.MolWt(m), 1))
    if "cLogP" in extra: col("cLogP", lambda m: round(Descriptors.MolLogP(m), 2))
    if "TPSA" in extra: col("TPSA", lambda m: round(Descriptors.TPSA(m), 1))
    if "HBD" in extra: col("HBD", Lipinski.NumHDonors)
    if "HBA" in extra: col("HBA", Lipinski.NumHAcceptors)
    if "Rot. bonds" in extra: col("Rot. bonds", Lipinski.NumRotatableBonds)
    if "Lipinski violations" in extra:
        col("Lipinski violations", lambda m: int(Descriptors.MolWt(m) > 500) + int(Descriptors.MolLogP(m) > 5) + int(Lipinski.NumHDonors(m) > 5) + int(Lipinski.NumHAcceptors(m) > 10))
    if "Scaffold" in extra: col("Scaffold", lambda m: Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) or "acyclic")
    if "Similarity to selected" in extra and ref is not None:
        out["Similarity to selected"] = [round(DataStructs.TanimotoSimilarity(ref, fps[k]), 2) for k in out["_ligand_key"]]
    return out


def structure_grid(R: dict, tab: pd.DataFrame, per_row: int = 3, size=(230, 170)) -> bytes:
    """RDKit Draw.MolsToGridImage of the listed compounds: legend = rank, id, CNN pose score, warhead; warhead atoms highlighted."""
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    prep = R["prepared"].set_index("ligand_key", drop=False)
    mols, legends, hl = [], [], []
    for _, r in tab.iterrows():
        p = prep.loc[r["_ligand_key"]]; m = Chem.RemoveHs(Chem.SDMolSupplier(p.ligand_sdf, removeHs=False)[0]); AllChem.Compute2DCoords(m)
        hits = hits_of(p); mols.append(m); hl.append(sorted({a for h in hits for a in h["match_atoms"]}))
        legends.append(f"{int(r['Rank'])}. {r['Name']}  |  CNN {r['CNN pose']:.2f}" + (f"  |  {r['Warhead']}" if r["Warhead"] != "—" else ""))
    return Draw.MolsToGridImage(mols, molsPerRow=per_row, subImgSize=size, legends=legends, highlightAtomLists=hl, returnPNG=True)
