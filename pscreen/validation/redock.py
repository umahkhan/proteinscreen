"""Redocking / pose-recovery benchmark on public RCSB complexes.

Ligand input is built from the CCD SMILES with a fresh RDKit conformer (NOT the crystal coordinates). The search box is
centred on the crystal ligand. RMSD is heavy-atom, symmetry-aware, no re-alignment. Success = RMSD <= 2.0 A.
Reported separately for the top pose under each GNINA signal, plus the best of all retained poses (an oracle, not a method).
LIMITATION: GNINA's CNN scorers were trained on PDBbind/CrossDocked-derived data; classic complexes below are very likely in that
training data, so success rates are optimistic and are NOT an estimate of performance on novel targets.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
from ..docking import DockConfig, run_screen, Store
from ..prep import prepare_library
from ..pockets import from_known_ligand
from .structs import fetch_pdb, ccd_smiles, split_complex, ref_ligand_mol, symmetric_rmsd, near_chains

NOT_LIGANDS = {"HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "CL", "NA", "K", "MG", "ZN", "CA", "MN", "FE", "NI", "CU", "CD", "IOD", "BR", "DMS",
               "FMT", "TRS", "MPD", "PGE", "PG4", "BME", "IMD", "NO3", "SCN", "EPE", "MES", "CIT", "TAR", "ACY", "NH4", "HG", "CO"}


def detect_ligand(pdb_path: Path, min_heavy: int = 10) -> tuple[str, str, int]:
    groups: dict[tuple, int] = {}
    for l in open(pdb_path):
        if l.startswith("HETATM") and l[17:20].strip() not in NOT_LIGANDS and l[76:78].strip() != "H":
            groups[(l[17:20].strip(), l[21], int(l[22:26]))] = groups.get((l[17:20].strip(), l[21], int(l[22:26])), 0) + 1
    cand = [(n, k) for k, n in groups.items() if n >= min_heavy]
    if not cand: raise ValueError(f"{pdb_path.stem}: no ligand with >= {min_heavy} heavy atoms")
    n, k = max(cand)
    return k


def build_task(pdb_id: str, workdir: Path, lig_resname: str | None = None, chains: list[str] | None = None) -> dict:
    pdb = fetch_pdb(pdb_id)
    name, chain, num = detect_ligand(pdb) if lig_resname is None else (lig_resname, None, None)
    chains = chains or near_chains(pdb, name, chain, num)
    parts = split_complex(pdb, chains, name, chain, num, workdir / pdb_id.upper())
    return {"pdb": pdb_id.upper(), "lig": name, "chains": "".join(chains), "smiles": ccd_smiles(name), **parts}


def run_redock(pdb_ids: list[str], outdir: str | Path, cfg: DockConfig, overrides: dict | None = None, progress=None) -> pd.DataFrame:
    outdir = Path(outdir).resolve(); outdir.mkdir(parents=True, exist_ok=True)
    overrides = overrides or {}
    tasks, pockets, lib_rows = [], [], []
    for pid in pdb_ids:
        try:
            t = build_task(pid, outdir / "inputs", **overrides.get(pid, {}))
        except Exception as e:
            tasks.append({"pdb": pid, "error": f"input build failed: {e}"}); continue
        pk = from_known_ligand(f"{t['pdb']}_site", str(t["receptor"]), str(t["ligand_ref"]), padding=6.0, min_size=20.0)
        pk.receptor_path = str(t["receptor"])
        pockets.append(pk); lib_rows.append({"id": t["pdb"], "smiles": t["smiles"]}); tasks.append(t)
    pd.DataFrame(lib_rows).to_csv(outdir / "ligands.csv", index=False)
    prepared = prepare_library(outdir / "ligands.csv", outdir)
    pairs = {(t["pdb"], f"{t['pdb']}_site") for t in tasks if "error" not in t}
    store = run_screen(prepared, pockets[0].receptor_path, pockets, cfg, outdir, progress=progress, pairs=pairs)
    return score_redock(outdir, tasks, prepared, store)


def score_redock(outdir: Path, tasks: list[dict], prepared: pd.DataFrame, store: Store) -> pd.DataFrame:
    poses = store.poses(); jobs = store.jobs().set_index("ligand_key")
    rows = []
    for t in tasks:
        if "error" in t:
            rows.append({"pdb": t["pdb"], "status": "input_failed", "reason": t["error"]}); continue
        row = {"pdb": t["pdb"], "ligand": t["lig"], "heavy_atoms": None, "status": "ok", "reason": None}
        prow = prepared[prepared.ligand_key == t["pdb"]]
        if prow.empty or prow.iloc[0].status != "ok":
            row.update(status="prep_failed", reason=None if prow.empty else prow.iloc[0].reason); rows.append(row); continue
        row["heavy_atoms"] = int(prow.iloc[0].heavy_atoms)
        g = poses[poses.ligand_key == t["pdb"]]
        if g.empty:
            row.update(status="docking_failed", reason=jobs.loc[t["pdb"], "error"] if t["pdb"] in jobs.index else "not run"); rows.append(row); continue
        ref = ref_ligand_mol(t["ligand_ref"])
        sup = [m for m in Chem.SDMolSupplier(g.iloc[0].pose_sdf, removeHs=False, sanitize=False)]
        g = g.assign(rmsd=[symmetric_rmsd(ref, Chem.RemoveHs(sup[int(i)], sanitize=False)) for i in g.pose_idx])
        if g.rmsd.isna().all():
            row.update(status="rmsd_unmapped", reason="docked ligand graph does not match crystal ligand (CCD SMILES vs deposited atoms)"); rows.append(row); continue
        row["n_poses"] = len(g)
        for name, col, asc in [("top1_vina_like", "minimized_affinity", True), ("top1_cnn_pose", "cnn_score", False),
                               ("top1_cnn_vs", "cnn_vs", False), ("top1_cnn_affinity", "cnn_affinity", False)]:
            r = g.sort_values(col, ascending=asc).iloc[0]; row[f"rmsd_{name}"] = round(float(r.rmsd), 2)
        row["rmsd_best_any_pose_oracle"] = round(float(g.rmsd.min()), 2)
        rows.append(row)
    res = pd.DataFrame(rows); res.to_csv(outdir / "redock_results.csv", index=False)
    return res


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0: return (float("nan"), float("nan"))
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return round((c - h) / d, 2), round((c + h) / d, 2)


def summarize(res: pd.DataFrame, cutoff: float = 2.0) -> dict:
    ok = res[res.status == "ok"]
    out = {"n_targets": len(res), "n_scored": len(ok), "n_unscored": len(res) - len(ok), "cutoff_A": cutoff}
    for c in [c for c in ok.columns if c.startswith("rmsd_")]:
        k = int((ok[c] <= cutoff).sum())
        out[c.replace("rmsd_", "success_")] = f"{k}/{len(ok)} = {k/len(ok):.2f} (Wilson 95% CI {wilson(k, len(ok))})" if len(ok) else None
    return out
