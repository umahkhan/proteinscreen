"""Compound library ingestion, standardization, dedup, descriptors, 3D prep, warhead annotation.

Every input record yields exactly one output row. Nothing is dropped silently: invalid and duplicate
records stay in the table with a status and reason, and keep their original identifier.
"""
from __future__ import annotations
import json, re
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from .warheads import find_warheads

RDLogger.DisableLog("rdApp.*")  # we collect our own warnings

ID_COLS = ("id", "compound_id", "name", "mol_id", "molecule_id", "title")
SMILES_COLS = ("smiles", "canonical_smiles", "smi", "isomeric_smiles")


@dataclass
class RawRecord:
    original_id: str
    source_index: int
    mol: Chem.Mol | None
    input_smiles: str | None
    error: str | None = None


def read_library(path: str | Path) -> list[RawRecord]:
    """Read .smi/.txt/.csv/.tsv (SMILES) or .sdf. Unparseable entries become records with `error` set."""
    path = Path(path)
    ext = path.suffix.lower()
    recs: list[RawRecord] = []
    if ext in (".sdf", ".sd"):
        suppl = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        for i, m in enumerate(suppl):
            if m is None:
                recs.append(RawRecord(f"record_{i+1}", i, None, None, "SDF record could not be parsed"))
                continue
            oid = next((m.GetProp(p) for p in ("ID", "id", "Name", "name", "_Name") if m.HasProp(p) and m.GetProp(p).strip()), f"record_{i+1}")
            recs.append(RawRecord(oid.strip(), i, m, None))
    elif ext in (".csv", ".tsv"):
        df = pd.read_csv(path, sep="\t" if ext == ".tsv" else ",", dtype=str, keep_default_na=False)
        low = {c.lower(): c for c in df.columns}
        sc = next((low[c] for c in SMILES_COLS if c in low), None)
        if sc is None:
            raise ValueError(f"No SMILES column found (looked for {SMILES_COLS}); columns: {list(df.columns)}")
        ic = next((low[c] for c in ID_COLS if c in low), None)
        for i, row in df.iterrows():
            _add_smiles(recs, row[ic].strip() if ic and row[ic].strip() else f"record_{i+1}", i, row[sc].strip())
    elif ext in (".smi", ".txt", ".smiles"):
        i = 0
        for line in path.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            _add_smiles(recs, parts[1].strip() if len(parts) > 1 else f"record_{i+1}", i, parts[0])
            i += 1
    else:
        raise ValueError(f"Unsupported library format: {ext}")
    return recs


def _add_smiles(recs, oid, idx, smi):
    if not smi:
        recs.append(RawRecord(oid, idx, None, smi, "Empty SMILES"))
        return
    m = Chem.MolFromSmiles(smi)
    recs.append(RawRecord(oid, idx, m, smi, None if m is not None else "SMILES could not be parsed/sanitized by RDKit"))


def _standardize(mol: Chem.Mol) -> tuple[Chem.Mol, list[str]]:
    warns = []
    m = rdMolStandardize.Cleanup(mol)
    if len(Chem.GetMolFrags(m)) > 1:
        m = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(m)
        warns.append("multiple fragments: kept largest (counter-ions/solvent removed)")
    before = Chem.GetFormalCharge(m)
    m2 = rdMolStandardize.Uncharger().uncharge(m)
    if Chem.GetFormalCharge(m2) != before:
        warns.append("neutralized charges (protonation state at assay pH not modeled)")
    m = m2
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return m, warns


def descriptors(mol: Chem.Mol) -> dict:
    return {
        "mol_weight": round(Descriptors.MolWt(mol), 3),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "clogp": round(Crippen.MolLogP(mol), 3),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rot_bonds": Lipinski.NumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "n_stereocenters_undefined": sum(1 for _, t in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False) if t == "?"),
    }


def embed_3d(mol: Chem.Mol, seed: int = 42, max_attempts: int = 10) -> tuple[Chem.Mol | None, str | None]:
    """Single ETKDGv3 + MMFF conformer with explicit H. GNINA does flexible docking, so this is a starting geometry only."""
    mh = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed; ps.maxIterations = max_attempts * 100
    if AllChem.EmbedMolecule(mh, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mh, ps) != 0:
            return None, "3D embedding failed"
    note = None
    try:
        if AllChem.MMFFOptimizeMolecule(mh, maxIters=2000) == 1:
            note = "MMFF optimization did not fully converge"
    except Exception:
        note = "MMFF parameters unavailable; unoptimized ETKDG geometry"
    return mh, note


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)[:80]


def prepare_library(path: str | Path, outdir: str | Path, write_3d: bool = True) -> pd.DataFrame:
    """Full preparation. Returns one row per input record; writes prepared.csv and ligands/*.sdf under outdir."""
    outdir = Path(outdir); (outdir / "ligands").mkdir(parents=True, exist_ok=True)
    recs = read_library(path)
    rows, seen = [], {}
    used_ids: dict[str, int] = {}
    for r in recs:
        # keep original ID verbatim; a unique ligand_key is derived for filenames
        n = used_ids.get(r.original_id, 0); used_ids[r.original_id] = n + 1
        key = _safe(r.original_id) + (f"__{n}" if n else "")
        row = {"compound_id": r.original_id, "ligand_key": key, "source_index": r.source_index,
               "input_smiles": r.input_smiles, "status": "ok", "warnings": [], "reason": None}
        if r.mol is None:
            row.update(status="invalid", reason=r.error); rows.append(row); continue
        try:
            mol = r.mol
            if r.input_smiles is None:   # SDF path: sanitize now, report failure
                Chem.SanitizeMol(mol)
            std, w = _standardize(mol)
            row["warnings"] += w
            ik = Chem.MolToInchiKey(std)
            row["canonical_smiles"] = Chem.MolToSmiles(std)
            row["inchikey"] = ik
            if r.input_smiles is None:
                row["input_smiles"] = Chem.MolToSmiles(mol)
        except Exception as e:
            row.update(status="invalid", reason=f"standardization failed: {e}"); rows.append(row); continue
        if not ik:
            row.update(status="invalid", reason="could not compute InChIKey"); rows.append(row); continue
        if ik in seen:
            row.update(status="duplicate", reason=f"duplicate of {seen[ik]}", duplicate_of=seen[ik]); rows.append(row); continue
        seen[ik] = r.original_id
        row.update(descriptors(std))
        hits = find_warheads(std)
        row["warheads"] = json.dumps(hits)
        row["warhead_classes"] = ";".join(sorted({h["warhead"] for h in hits}))
        row["reactive_atom_indices"] = json.dumps({h["warhead"] + f"#{i}": h["reactive_atoms"] for i, h in enumerate(hits)})
        if not hits:
            pass  # remains eligible for noncovalent ranking
        if std.GetNumHeavyAtoms() > 80:
            row["warnings"].append("very large molecule (>80 heavy atoms); docking may be slow/unreliable")
        if write_3d:
            m3, note = embed_3d(std)
            if m3 is None:
                row.update(status="prep_failed", reason=note)
            else:
                if note: row["warnings"].append(note)
                m3.SetProp("_Name", r.original_id); m3.SetProp("ligand_key", key)
                p = outdir / "ligands" / f"{key}.sdf"
                with Chem.SDWriter(str(p)) as wtr:
                    wtr.write(m3)
                row["ligand_sdf"] = str(p)
        rows.append(row)
    df = pd.DataFrame(rows)
    df["warnings"] = df["warnings"].apply(lambda w: "; ".join(w))
    df.to_csv(outdir / "prepared.csv", index=False)
    return df


def summarize(df: pd.DataFrame) -> dict:
    return {"total": len(df), **df["status"].value_counts().to_dict(),
            "with_warhead": int(df["warhead_classes"].fillna("").ne("").sum()) if "warhead_classes" in df else 0}
