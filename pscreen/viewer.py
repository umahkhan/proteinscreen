"""3Dmol.js scene builder (via py3Dmol). All coordinates come from real files: receptor PDB, GNINA pose SDFs,
experimental reference ligand. Nothing is generated or illustrative.

Colour code:  receptor = grey cartoon | pocket residues = orange sticks | PREDICTED poses = cyan/magenta/... carbons
              EXPERIMENTAL reference ligand = green, thin sticks + 'EXPERIMENTAL' label | warhead atoms = red spheres
              reactive atom = large red sphere | target residue = yellow sticks | distance = dashed yellow line
"""
from __future__ import annotations
from rdkit import Chem

POSE_COLORS = ["cyan", "magenta", "orange", "lightblue", "pink", "wheat", "violet", "salmon"]


def pose_block(pose_sdf: str, pose_idx: int) -> Chem.Mol:
    return [m for m in Chem.SDMolSupplier(pose_sdf, removeHs=False, sanitize=False)][pose_idx]
