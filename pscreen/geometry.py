"""Pose -> ligand atom mapping and reactive-atom / target-residue geometry.

Heuristic thresholds (assumptions, not validated cutoffs):
  <= 4.0 A  'plausible'   reactive atom within contact range of the nucleophile in a NONcovalent pose
  <= 6.0 A  'marginal'
  >  6.0 A  'unlikely'
These describe pose geometry only; they say nothing about reactivity, kinetics, potency or selectivity.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from rdkit import Chem
from Bio.PDB import MMCIFParser, PDBParser

# residue -> candidate nucleophilic atoms (nearest is used)
NUCLEOPHILE_ATOMS = {"CYS": ["SG"], "LYS": ["NZ"], "SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"],
                     "HIS": ["NE2", "ND1"], "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
PLAUSIBLE_A, MARGINAL_A = 4.0, 6.0


def _skeleton(mol: Chem.Mol) -> tuple[Chem.Mol, list[int]]:
    """Element-only, single-bond heavy-atom graph (robust to bond-order/aromaticity differences). Returns (skel, orig_idx)."""
    rw = Chem.RWMol(); idx = {}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            idx[a.GetIdx()] = rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in idx and j in idx:
            rw.AddBond(idx[i], idx[j], Chem.BondType.SINGLE)
    return rw.GetMol(), list(idx)


def map_ligand_to_pose(ligand: Chem.Mol, pose: Chem.Mol) -> dict[int, int] | None:
    """Map ligand heavy-atom indices (prepared SDF order) -> pose atom indices. None if graphs don't match."""
    ls, l_orig = _skeleton(ligand); ps, p_orig = _skeleton(pose)
    if ls.GetNumAtoms() != ps.GetNumAtoms():
        return None
    m = ps.GetSubstructMatch(ls)
    if not m or len(m) != ls.GetNumAtoms():
        return None
    return {l_orig[qi]: p_orig[ti] for qi, ti in enumerate(m)}


def load_structure(path: str):
    p = MMCIFParser(QUIET=True) if str(path).lower().endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    return p.get_structure("rec", path)


def parse_residue(spec: str) -> tuple[str, int, str]:
    """'A:CYS145' or 'A:145' -> (chain, resnum, resname or '')"""
    ch, rest = spec.split(":")
    name = "".join(c for c in rest if c.isalpha()).upper()
    num = int("".join(c for c in rest if c.isdigit() or c == "-"))
    return ch, num, name


def residue_atoms(structure, spec: str) -> dict[str, tuple[np.ndarray, str]]:
    """Return {atom_name: (xyz, resname)} for the target residue."""
    ch, num, name = parse_residue(spec)
    for model in structure:
        if ch in model:
            for res in model[ch]:
                if res.id[1] == num and res.id[0].strip() == "":
                    if name and res.get_resname() != name:
                        raise ValueError(f"{spec}: residue {num} on chain {ch} is {res.get_resname()}, not {name}")
                    return {a.get_id(): (a.coord.copy(), res.get_resname()) for a in res}
        break
    raise ValueError(f"Residue {spec} not found")


@dataclass
class Geometry:
    distance: float | None
    nucleophile_atom: str | None
    reactive_atom_pose_index: int | None
    label: str                    # plausible|marginal|unlikely|unmapped|not_applicable
    target_xyz: tuple | None = None
    reactive_xyz: tuple | None = None


def reactive_geometry(ligand: Chem.Mol, pose: Chem.Mol, reactive_atoms: list[int], target_atoms: dict) -> Geometry:
    """Min distance over (reactive atom, candidate nucleophile atom) pairs."""
    resname = next(iter(target_atoms.values()))[1]
    cands = [n for n in NUCLEOPHILE_ATOMS.get(resname, []) if n in target_atoms]
    if not cands:
        return Geometry(None, None, None, "not_applicable")
    mp = map_ligand_to_pose(ligand, pose)
    if mp is None or any(a not in mp for a in reactive_atoms):
        return Geometry(None, None, None, "unmapped")
    xyz = pose.GetConformer().GetPositions()
    best = None
    for ra in reactive_atoms:
        for n in cands:
            d = float(np.linalg.norm(xyz[mp[ra]] - target_atoms[n][0]))
            if best is None or d < best[0]:
                best = (d, n, mp[ra], tuple(xyz[mp[ra]]), tuple(target_atoms[n][0]))
    d, n, pi, rxyz, txyz = best
    label = "plausible" if d <= PLAUSIBLE_A else "marginal" if d <= MARGINAL_A else "unlikely"
    return Geometry(round(d, 2), n, int(pi), label, tuple(float(x) for x in txyz), tuple(float(x) for x in rxyz))
