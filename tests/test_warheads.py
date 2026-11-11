import pytest
from rdkit import Chem
from pscreen.warheads import WARHEADS, find_warheads


def names(smi):
    m = Chem.MolFromSmiles(smi)
    assert m is not None, smi
    return {h["warhead"] for h in find_warheads(m)}


@pytest.mark.parametrize("w", WARHEADS, ids=lambda w: w.name)
def test_positives_match(w):
    assert w.positives, f"{w.name} has no positive controls"
    for smi in w.positives:
        assert w.name in names(smi), f"{w.name} should match {smi}"


@pytest.mark.parametrize("w", WARHEADS, ids=lambda w: w.name)
def test_negatives_do_not_match(w):
    assert w.negatives, f"{w.name} has no negative controls"
    for smi in w.negatives:
        assert w.name not in names(smi), f"{w.name} should NOT match {smi}"


def test_reactive_atom_is_beta_carbon_of_acrylamide():
    m = Chem.MolFromSmiles("C=CC(=O)Nc1ccccc1")
    (h,) = [h for h in find_warheads(m) if h["warhead"] == "acrylamide"]
    assert h["reactive_atoms"] == [0]
    assert m.GetAtomWithIdx(0).GetSymbol() == "C" and m.GetAtomWithIdx(0).GetTotalNumHs() == 2


def test_no_warhead_compound_returns_empty():
    assert find_warheads(Chem.MolFromSmiles("c1ccccc1CCO")) == []
