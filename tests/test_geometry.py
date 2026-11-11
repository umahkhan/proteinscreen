import json
import numpy as np
from rdkit import Chem
from pscreen.geometry import map_ligand_to_pose
from pscreen.warheads import find_warheads


def test_mapping_survives_atom_reordering_and_dropped_hydrogens():
    lig = Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1")
    hit = [h for h in find_warheads(lig) if h["warhead"] == "acrylamide"][0]
    # shuffle atom order (as GNINA/OpenBabel does) and use a bond-order-free, unsanitized pose
    order = list(range(lig.GetNumAtoms()))[::-1]
    pose = Chem.RenumberAtoms(lig, order)
    mp = map_ligand_to_pose(lig, pose)
    assert mp is not None
    ra = hit["reactive_atoms"][0]
    assert pose.GetAtomWithIdx(mp[ra]).GetTotalNumHs() == 2 and pose.GetAtomWithIdx(mp[ra]).GetSymbol() == "C"


def test_mismatched_graph_returns_none():
    assert map_ligand_to_pose(Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("CCN")) is None
