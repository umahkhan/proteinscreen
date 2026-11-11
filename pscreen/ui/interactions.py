"""Geometric protein-ligand contact analysis of a docked pose (heuristics only - not an energy calculation).

Types: 'H-bond candidate' (N/O...N/O <= 3.5 A), 'ionic/polar' (charged side-chain atom to ligand N/O <= 4.0 A),
'hydrophobic' (C...C <= 4.0 A), 'other close contact' (<= 4.0 A). Clash flag: any non-covalent heavy-atom pair < CLASH_A.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from rdkit import Chem
from scipy.spatial import cKDTree
from ..geometry import load_structure

CLASH_A = 2.5
CHARGED = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("HIS", "ND1"), ("HIS", "NE2")}


def _atoms(receptor: str):
    rows = []
    for model in load_structure(receptor):
        for ch in model:
            for r in ch:
                if r.id[0].strip() == "":
                    for a in r:
                        if a.element != "H": rows.append((ch.id, r.id[1], r.get_resname(), a.get_id(), a.element, a.coord.copy()))
        break
    return rows


_cache: dict = {}


def analyse(receptor: str, pose: Chem.Mol, cutoff: float = 4.0, exclude_pairs: set | None = None) -> tuple[pd.DataFrame, dict]:
    """Returns (per-residue table, summary dict). exclude_pairs: {(pose_atom_idx, 'CHAIN:RESNUM:ATOM')} to ignore (modelled covalent bond)."""
    if receptor not in _cache: _cache[receptor] = _atoms(receptor)
    P = _cache[receptor]; pxyz = np.array([a[5] for a in P])
    lig = [a for a in pose.GetAtoms() if a.GetAtomicNum() > 1]
    lxyz = pose.GetConformer().GetPositions()[[a.GetIdx() for a in lig]]
    tree = cKDTree(pxyz); rows, clashes, min_d = {}, 0, 99.0
    for li, a in enumerate(lig):
        for pi in tree.query_ball_point(lxyz[li], cutoff):
            ch, num, rn, an, el, xyz = P[pi]; d = float(np.linalg.norm(lxyz[li] - xyz))
            if exclude_pairs and (a.GetIdx(), f"{ch}:{num}:{an}") in exclude_pairs: continue
            min_d = min(min_d, d)
            if d < CLASH_A: clashes += 1
            lel = a.GetSymbol()
            if lel in "NO" and el in "NO" and d <= 3.5: t = "H-bond candidate"
            elif lel in "NO" and (rn, an) in CHARGED and d <= 4.0: t = "ionic/polar"
            elif lel == "C" and el == "C" and d <= 4.0: t = "hydrophobic"
            else: t = "other close contact"
            k = (ch, num, rn); r = rows.setdefault(k, {"residue": f"{ch}:{rn}{num}", "types": set(), "closest_A": 99.0, "n_atom_pairs": 0, "closest_pair": ""})
            r["types"].add(t); r["n_atom_pairs"] += 1
            if d < r["closest_A"]: r["closest_A"], r["closest_pair"] = round(d, 2), f"{lel}{a.GetIdx()}…{an}"
    df = pd.DataFrame([{**r, "types": ", ".join(sorted(r["types"]))} for r in rows.values()])
    if len(df): df = df.sort_values("closest_A").reset_index(drop=True)
    summ = {"residues_in_contact": len(df), "hbond_candidate_residues": int(df.types.str.contains("H-bond").sum()) if len(df) else 0,
            "hydrophobic_residues": int(df.types.str.contains("hydrophobic").sum()) if len(df) else 0, "clashes_lt_2.5A": clashes,
            "closest_protein_ligand_A": round(min_d, 2) if min_d < 99 else None}
    return df, summ
