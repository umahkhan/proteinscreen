import pytest
from rdkit import Chem
from pscreen.covalent import covalent_candidates

def canon(m): return Chem.MolToSmiles(m)

# (input SMILES, warhead class, expected adduct SMILES (canonical), attach atom element, attach total H)
CASES = [
    ("C=CC(=O)Nc1ccccc1", "acrylamide", "CCC(=O)Nc1ccccc1", "C", 3),
    ("ClCC(=O)NC", "chloroacetamide", "CC(=O)NC", "C", 3),
    ("C#CC(=O)NC", "alkynamide_or_ynone", "C=CC(=O)NC", "C", 2),
    ("O=Cc1ccccc1", "aldehyde", "OCc1ccccc1", "C", 2),
    ("N#Cc1ccccc1", "nitrile", "N=Cc1ccccc1", "C", 1),
    ("O=S(=O)(F)c1ccccc1", "sulfonyl_fluoride", "O=[SH](=O)c1ccccc1", "S", 1),
    ("O=C1C=CC(=O)N1C", "maleimide", "CN1C(=O)CCC1=O", "C", 2),
    ("CN=C=S", "isothiocyanate", "CNC=S", "C", 1),
    ("Clc1ncccn1", "snar_2_halopyrimidine", "c1ccncn1", "C", 1),
]

@pytest.mark.parametrize("smi,wh,exp,el,nh", CASES)
def test_adduct_forms(smi, wh, exp, el, nh):
    forms, _ = covalent_candidates(Chem.MolFromSmiles(smi))
    f = next(f for f in forms if f.warhead == wh)
    assert canon(f.product) == canon(Chem.MolFromSmiles(exp)), canon(f.product)
    a = f.product.GetAtomWithIdx(f.attach_idx)
    assert a.GetSymbol() == el and a.GetTotalNumHs() == nh
    # the SMARTS must match the attachment atom in the product
    patt = Chem.MolFromSmarts(f.gnina_smarts)
    assert patt is not None and f.attach_idx in {m[0] for m in f.product.GetSubstructMatches(patt)}

def test_epoxide_gives_two_regioisomeric_forms():
    forms, _ = covalent_candidates(Chem.MolFromSmiles("CC1OC1c1ccccc1"))
    ep = [f for f in forms if f.warhead == "epoxide"]
    # attack at the methyl-bearing carbon -> PhCH(OH)CH2CH3; attack at the phenyl-bearing carbon -> PhCH2CH(OH)CH3
    assert {canon(f.product) for f in ep} == {canon(Chem.MolFromSmiles("CCC(O)c1ccccc1")), canon(Chem.MolFromSmiles("CC(O)Cc1ccccc1"))}
    assert len({f.attach_idx for f in ep}) == 2

def test_unsupported_classes_reported_not_dropped():
    forms, unsup = covalent_candidates(Chem.MolFromSmiles("CSSC"))
    assert forms == [] and unsup == ["disulfide"]

def test_no_warhead_no_forms():
    assert covalent_candidates(Chem.MolFromSmiles("c1ccccc1CCO")) == ([], [])

def test_source_to_product_mapping_after_atom_removal():
    mol = Chem.MolFromSmiles("ClCC(=O)NC")
    f = covalent_candidates(mol)[0][0]
    assert 0 not in f.source_to_product                          # Cl removed
    assert f.product.GetAtomWithIdx(f.source_to_product[3]).GetSymbol() == "O"
