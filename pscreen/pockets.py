"""Pocket definition. Known-pocket and user-box paths (no prediction). Predicted-pocket tools plug in via the same Pocket type.

Docking box method (documented, user-adjustable):
  center = centroid of the reference ligand (or user coordinates / centroid of selected residues' atoms)
  size   = ligand extent along each axis + 2*padding (default padding 6 A), floored at min_size (default 18 A)
Residue-level cryptic-pocket probabilities do NOT define docking coordinates by themselves; see pocket tools for their own box logic.
"""
from __future__ import annotations
import numpy as np
from Bio.PDB import NeighborSearch
from .docking import Pocket
from .geometry import load_structure


def _ligand_xyz_pdb(path: str) -> np.ndarray:
    xyz = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in open(path)
           if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"]
    if not xyz: raise ValueError(f"No heavy atoms read from {path}")
    return np.array(xyz)


def residues_near(receptor: str, xyz: np.ndarray, cutoff: float = 5.0) -> list[str]:
    st = load_structure(receptor)
    atoms = [a for a in st.get_atoms() if a.get_parent().id[0].strip() == ""]
    ns = NeighborSearch(atoms)
    found = set()
    for p in xyz:
        for a in ns.search(p, cutoff):
            r = a.get_parent(); found.add((r.get_parent().id, r.id[1], r.get_resname()))
    return [f"{c}:{n}{i}" for c, i, n in sorted(found)]


def from_known_ligand(pocket_id: str, receptor: str, ligand_pdb: str, padding: float = 6.0, min_size: float = 18.0) -> Pocket:
    xyz = _ligand_xyz_pdb(ligand_pdb)
    center = xyz.mean(0); ext = xyz.max(0) - xyz.min(0)
    size = np.maximum(ext + 2 * padding, min_size)
    return Pocket(pocket_id, tuple(map(float, center.round(3))), tuple(map(float, size.round(2))), "known_ligand",
                  receptor_path=receptor, residues=residues_near(receptor, xyz),
                  metrics={"reference_ligand": ligand_pdb, "evidence": "experimental co-crystal ligand position"})


def from_user_box(pocket_id: str, receptor: str, center, size=(22, 22, 22)) -> Pocket:
    c = np.array(center, float); h = np.array(size, float) / 2
    near = residues_near(receptor, c[None, :], cutoff=float(h.max()))
    return Pocket(pocket_id, tuple(c), tuple(size), "user_box", receptor_path=receptor, residues=near,
                  metrics={"evidence": "user-supplied coordinates"})
