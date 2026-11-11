import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pscreen.validation.structs import symmetric_rmsd


def _mol(smi, seed=1):
    m = Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m, randomSeed=seed); return Chem.RemoveHs(m)


def test_identity_is_zero_and_symmetric_swap_is_zero():
    ref = _mol("c1ccccc1C(=O)O")
    assert symmetric_rmsd(ref, Chem.Mol(ref)) < 1e-6
    # relabel atoms (graph-equivalent, different order) -> still 0
    order = list(range(ref.GetNumAtoms()))[::-1]
    assert symmetric_rmsd(ref, Chem.RenumberAtoms(ref, order)) < 1e-6


def test_carboxylate_oxygen_swap_handled():
    ref = _mol("CC(=O)O"); pose = Chem.Mol(ref)
    conf = pose.GetConformer(); p = [np.array(conf.GetAtomPosition(i)) for i in range(pose.GetNumAtoms())]
    # swap coordinates of the two oxygens (indices 2 and 3): symmetric RMSD must still be ~0
    conf.SetAtomPosition(2, p[3].tolist()); conf.SetAtomPosition(3, p[2].tolist())
    assert symmetric_rmsd(ref, pose) < 1e-6


def test_translation_gives_known_rmsd():
    ref = _mol("CCO"); pose = Chem.Mol(ref); conf = pose.GetConformer()
    for i in range(pose.GetNumAtoms()):
        x = conf.GetAtomPosition(i); conf.SetAtomPosition(i, (x.x + 3.0, x.y, x.z))
    assert abs(symmetric_rmsd(ref, pose) - 3.0) < 1e-6


def test_mismatch_returns_none():
    assert symmetric_rmsd(_mol("CCO"), _mol("CCN")) is None
