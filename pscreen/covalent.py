"""Bound-state ligand construction for GNINA covalent docking.

GNINA (1.3) expects the ligand already in its post-reaction chemistry (ring opened, leaving group removed, unsaturation
saturated) with the bonding atom specified by SMARTS; the receptor atom is given by chain:resnum:atom. We therefore
convert a recognised warhead (see warheads.py) into its adduct form using the convention
    "the bond to the protein is represented by an implicit H on the attachment atom"
(e.g. acrylamide C=CC(=O)N -> CCC(=O)N, attachment = the former beta carbon, now CH3).
This is a modelling convention for pose generation. It is NOT a reaction-energy calculation: it says nothing about whether
or how fast the reaction occurs, potency, or selectivity.

Classes without a defined transform (boronic acid, disulfide) are reported as unsupported, never silently skipped.
"""
from __future__ import annotations
from dataclasses import dataclass
from rdkit import Chem
from .warheads import find_warheads

MICHAEL = {"acrylamide", "acrylate_or_enone", "vinyl_sulfone_or_sulfonamide", "nitroalkene", "cyanoacrylamide", "maleimide"}
ALKYLATING = {"chloroacetamide", "alpha_haloketone"}
SUPPORTED = MICHAEL | ALKYLATING | {"alkynamide_or_ynone", "snar_2_halopyrimidine", "epoxide", "aziridine", "isothiocyanate",
                                    "sulfonyl_fluoride", "aldehyde", "nitrile"}


@dataclass
class BoundForm:
    warhead: str
    product: Chem.Mol            # sanitized adduct-form molecule (heavy atoms + implicit H)
    attach_idx: int              # attachment atom index in `product`
    source_to_product: dict[int, int]   # original heavy-atom index -> product index (removed atoms absent)
    gnina_smarts: str            # SMARTS for GNINA --covalent_lig_atom_pattern
    note: str = ""


def _rw(mol):
    m = Chem.RWMol(Chem.Mol(mol)); Chem.Kekulize(m, clearAromaticFlags=False); return m


def _finish(rw: Chem.RWMol, removed: set[int], attach: int, warhead: str, src: Chem.Mol, note="") -> BoundForm | None:
    keep = [i for i in range(rw.GetNumAtoms()) if i not in removed]
    mapping = {old: new for new, old in enumerate(keep)}
    for i in sorted(removed, reverse=True):
        rw.RemoveAtom(i)
    m = rw.GetMol()
    for a in m.GetAtoms():
        a.SetNoImplicit(False); a.SetNumExplicitHs(0)       # let valence rules supply the H that stands in for the protein bond
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    new_attach = mapping[attach]
    return BoundForm(warhead, m, new_attach, mapping, attach_smarts(m, new_attach), note)


def attach_smarts(mol: Chem.Mol, idx: int) -> str:
    """Recursive SMARTS for the attachment atom: element + the elements of its heavy neighbours (bond types/H counts left free,
    so it survives re-protonation by GNINA's preprocessing). GNINA evaluates every matching atom."""
    a = mol.GetAtomWithIdx(idx)
    nbrs = "".join(f"(~[#{n.GetAtomicNum()}])" for n in a.GetNeighbors())
    return f"[$([#{a.GetAtomicNum()}]{nbrs})]"


def bound_forms(mol: Chem.Mol, hit: dict) -> list[BoundForm]:
    """All adduct forms for one warhead hit (epoxides/aziridines yield one per ring carbon). `mol` has no explicit H."""
    w, out = hit["warhead"], []
    if w not in SUPPORTED:
        return out
    for ra in hit["reactive_atoms"]:
        rw = _rw(mol); a = rw.GetAtomWithIdx(ra); removed: set[int] = set()
        if w in MICHAEL or w == "alkynamide_or_ynone":
            want = Chem.BondType.TRIPLE if w == "alkynamide_or_ynone" else Chem.BondType.DOUBLE
            b = next((b for b in a.GetBonds() if b.GetBondType() == want and b.GetOtherAtom(a).GetAtomicNum() == 6), None)
            if b is None: continue
            b.SetBondType(Chem.BondType.DOUBLE if want == Chem.BondType.TRIPLE else Chem.BondType.SINGLE)
        elif w in ALKYLATING or w in ("snar_2_halopyrimidine", "sulfonyl_fluoride"):
            hal = [n.GetIdx() for n in a.GetNeighbors() if n.GetSymbol() in (("F",) if w == "sulfonyl_fluoride" else ("F", "Cl", "Br", "I"))]
            if not hal: continue
            removed.add(hal[0])
        elif w in ("epoxide", "aziridine"):
            het = [n for n in a.GetNeighbors() if n.GetSymbol() in ("O", "N") and n.IsInRing()]
            if not het: continue
            rw.RemoveBond(ra, het[0].GetIdx())
        elif w == "isothiocyanate":
            b = next((b for b in a.GetBonds() if b.GetOtherAtom(a).GetSymbol() == "N" and b.GetBondType() == Chem.BondType.DOUBLE), None)
            if b is None: continue
            b.SetBondType(Chem.BondType.SINGLE)
        elif w == "aldehyde":
            b = next((b for b in a.GetBonds() if b.GetOtherAtom(a).GetSymbol() == "O" and b.GetBondType() == Chem.BondType.DOUBLE), None)
            if b is None: continue
            b.SetBondType(Chem.BondType.SINGLE)
        elif w == "nitrile":
            b = next((b for b in a.GetBonds() if b.GetBondType() == Chem.BondType.TRIPLE), None)
            if b is None: continue
            b.SetBondType(Chem.BondType.DOUBLE)
        bf = _finish(rw, removed, ra, w, mol)
        if bf: out.append(bf)
    return out


def covalent_candidates(mol: Chem.Mol) -> tuple[list[BoundForm], list[str]]:
    """Adduct forms for every recognised warhead, plus the list of recognised-but-unsupported classes."""
    forms, unsupported = [], []
    for h in find_warheads(mol):
        f = bound_forms(mol, h)
        if f: forms += f
        elif h["warhead"] not in SUPPORTED: unsupported.append(h["warhead"])
    return forms, sorted(set(unsupported))
