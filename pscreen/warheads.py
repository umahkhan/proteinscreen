"""Curated electrophilic warhead SMARTS library.

A match is a *chemical annotation only*. It does not establish reactivity, covalent engagement,
potency or selectivity. Reactive atoms are marked with atom-map numbers in the SMARTS (":1", ":2");
map numbers do not constrain matching, they only label which matched atoms are electrophilic centres.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from rdkit import Chem


@dataclass(frozen=True)
class Warhead:
    name: str
    smarts: str
    mechanism: str
    typical_residues: tuple[str, ...]
    reversible: bool
    note: str = ""
    positives: tuple[str, ...] = ()   # SMILES that must match
    negatives: tuple[str, ...] = ()   # SMILES that must not match
    _query: Chem.Mol = field(init=False, repr=False, compare=False, default=None)

    def __post_init__(self):
        q = Chem.MolFromSmarts(self.smarts)
        if q is None:
            raise ValueError(f"Invalid SMARTS for {self.name}: {self.smarts}")
        if not any(a.GetAtomMapNum() for a in q.GetAtoms()):
            raise ValueError(f"{self.name}: no reactive atom map number")
        object.__setattr__(self, "_query", q)

    @property
    def query(self) -> Chem.Mol:
        return self._query


WARHEADS: tuple[Warhead, ...] = (
    Warhead("acrylamide", "[CX3;!R;H1,H2:1]=[CX3]-C(=O)-[NX3]", "Michael addition (beta carbon)", ("Cys",), False,
            "Reactive atom is the beta carbon. Beta/alpha substitution strongly modulates reactivity.",
            ("C=CC(=O)Nc1ccccc1", "CC=CC(=O)NC", "C=C(C)C(=O)N(C)C"),
            ("CC(=O)NC", "O=C1C=CC(=O)N1", "O=C1NC=CC=C1")),
    Warhead("acrylate_or_enone", "[CX3;!R;H1,H2:1]=[CX3]-C(=O)-[O,#6]", "Michael addition (beta carbon)", ("Cys",), False,
            "Acrylate esters and vinyl/enone ketones. Typically more reactive than acrylamides.",
            ("C=CC(=O)OC", "C=CC(=O)c1ccccc1", "CC=CC(=O)C"),
            ("CC(=O)OC", "c1ccccc1C(=O)C")),
    Warhead("vinyl_sulfone_or_sulfonamide", "[CX3;!R;H1,H2:1]=[CX3]-[SX4](=O)(=O)", "Michael addition (beta carbon)", ("Cys", "Lys"), False,
            "",
            ("C=CS(=O)(=O)N", "C=CS(=O)(=O)c1ccccc1"),
            ("CS(=O)(=O)NC", "c1ccccc1S(=O)(=O)N")),
    Warhead("alkynamide_or_ynone", "[CX2:1]#[CX2]-C(=O)", "Michael addition (beta carbon)", ("Cys",), False,
            "Propiolamides/ynones.",
            ("C#CC(=O)N", "CC#CC(=O)Nc1ccccc1", "C#CC(=O)C"),
            ("C#CCN", "CC(=O)N")),
    Warhead("cyanoacrylamide", "[CX3:1]=[CX3](C#N)-C(=O)", "Reversible Michael addition", ("Cys",), True,
            "Reversible covalent warhead.",
            ("N#CC(=Cc1ccccc1)C(=O)N", "CC=C(C#N)C(=O)N"),
            ("N#CCC(=O)N", "C=CC(=O)N")),
    Warhead("nitroalkene", "[CX3:1]=[CX3]-[N+](=O)[O-]", "Michael addition", ("Cys",), True, "",
            ("C=C[N+](=O)[O-]", "c1ccccc1C=C[N+](=O)[O-]"),
            ("CC[N+](=O)[O-]", "c1ccccc1[N+](=O)[O-]")),
    Warhead("maleimide", "[CX3;H1:1]1=[CX3;H1:2]C(=O)[NX3]C1=O", "Michael addition", ("Cys",), False, "",
            ("O=C1C=CC(=O)N1", "O=C1C=CC(=O)N1C"),
            ("O=C1CCC(=O)N1", "C=CC(=O)N")),
    Warhead("chloroacetamide", "[CH2:1]([Cl,Br,I])-C(=O)-[NX3]", "SN2 alkylation", ("Cys",), False,
            "Reactive atom is the alpha carbon (leaving group departs).",
            ("ClCC(=O)Nc1ccccc1", "BrCC(=O)N(C)C"),
            ("CC(=O)NC", "ClC(Cl)C(=O)N", "ClCCN")),
    Warhead("alpha_haloketone", "[CH2:1]([Cl,Br,I])-C(=O)-[#6]", "SN2 alkylation", ("Cys", "His"), False, "",
            ("ClCC(=O)c1ccccc1", "BrCC(=O)C"),
            ("CC(=O)c1ccccc1", "ClCC(=O)N")),
    Warhead("epoxide", "[CX4:1]1O[CX4:2]1", "Ring opening (SN2)", ("Cys", "Asp", "Glu", "Thr"), False, "",
            ("C1CO1", "CC1OC1c1ccccc1"),
            ("C1CCOC1", "COC")),
    Warhead("aziridine", "[CX4:1]1[NX3][CX4:2]1", "Ring opening (SN2)", ("Cys",), False, "",
            ("C1CN1", "CC1CN1C(=O)C"),
            ("C1CCNC1", "CNC")),
    Warhead("isothiocyanate", "[NX2]=[CX2:1]=S", "Addition to nucleophile", ("Cys", "Lys"), True, "",
            ("CN=C=S", "S=C=Nc1ccccc1"),
            ("CNC(=S)N", "CSC#N")),
    Warhead("sulfonyl_fluoride", "[SX4:1](=O)(=O)F", "SuFEx", ("Lys", "Tyr", "Ser", "His"), False,
            "Includes fluorosulfates. Targets non-Cys residues.",
            ("O=S(=O)(F)c1ccccc1", "CS(=O)(=O)F", "O=S(=O)(F)Oc1ccccc1"),
            ("CS(=O)(=O)N", "CS(=O)(=O)Cl")),
    Warhead("aldehyde", "[CX3H1:1](=O)[#6]", "Reversible hemithioacetal / imine", ("Cys", "Lys", "Ser"), True, "",
            ("CC=O", "O=Cc1ccccc1"),
            ("CC(=O)C", "O=CN", "OC=O")),
    Warhead("nitrile", "[CX2:1](#[NX1])[#6]", "Reversible thioimidate formation", ("Cys", "Ser"), True,
            "Very common, usually non-reactive. Weak annotation only.",
            ("CC#N", "N#Cc1ccccc1"),
            ("CN", "CNC=O")),
    Warhead("boronic_acid", "[#6][BX3:1]([OX2])[OX2]", "Reversible covalent (Ser/Thr/Lys)", ("Ser", "Thr", "Lys"), True, "",
            ("OB(O)c1ccccc1", "CB(O)O", "CC1(C)OB(c2ccccc2)OC1(C)C"),
            ("CO", "c1ccccc1", "OB(O)O")),
    Warhead("disulfide", "[SX2:1][SX2:2]", "Thiol-disulfide exchange", ("Cys",), True, "",
            ("CSSC", "c1ccccc1SSc1ccccc1"),
            ("CSC", "CS")),
    Warhead("snar_2_halopyrimidine", "[c:1]([F,Cl])(:n):n", "SNAr", ("Cys",), False,
            "Low-confidence: reactivity depends strongly on ring substituents.",
            ("Clc1ncccn1", "Fc1nc(N)ccn1"),
            ("Clc1ccccc1", "Clc1ccncc1")),
)

_BY_NAME = {w.name: w for w in WARHEADS}


def get_warhead(name: str) -> Warhead:
    return _BY_NAME[name]


def find_warheads(mol: Chem.Mol) -> list[dict]:
    """Return all matches: {warhead, match_atoms, reactive_atoms, reversible, mechanism, ...}.

    Overlapping matches of different classes are all reported (e.g. cyanoacrylamide also matches nitrile
    and acrylamide). Identical (class, reactive atom set) duplicates from symmetric matches are collapsed.
    """
    out, seen = [], set()
    for w in WARHEADS:
        q = w.query
        maps = {a.GetIdx(): a.GetAtomMapNum() for a in q.GetAtoms() if a.GetAtomMapNum()}
        for match in mol.GetSubstructMatches(q, uniquify=True, maxMatches=1000):
            reactive = tuple(sorted(match[qi] for qi in maps))
            key = (w.name, reactive)
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "warhead": w.name, "match_atoms": list(match), "reactive_atoms": list(reactive),
                "mechanism": w.mechanism, "reversible": w.reversible,
                "typical_residues": list(w.typical_residues), "note": w.note,
            })
    return out
