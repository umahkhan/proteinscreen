"""2D structure depictions (RDKit) as PNG data URIs for table thumbnails and SVG for cards."""
from __future__ import annotations
import base64
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


# ---- depictions whose atom indices follow the prepared ligand SDF (== standardized-molecule order used for warhead detection)
_HL, _RX = (0.99, 0.80, 0.55), (0.92, 0.26, 0.20)      # warhead atoms: soft orange; reactive atom: red


def _ligand_mol(ligand_sdf: str) -> Chem.Mol:
    from rdkit.Chem import AllChem
    m = Chem.RemoveHs(Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]); AllChem.Compute2DCoords(m); return m


@lru_cache(maxsize=4096)
def compound_png(ligand_sdf: str, match_atoms: tuple = (), reactive_atoms: tuple = (), w: int = 170, h: int = 110) -> str:
    """Compound drawing with the warhead atoms shaded and the reactive atom in red."""
    m = _ligand_mol(ligand_sdf)
    d = rdMolDraw2D.MolDraw2DCairo(w, h); d.drawOptions().clearBackground = True
    cols = {i: _HL for i in match_atoms}; cols.update({i: _RX for i in reactive_atoms})
    rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=list(cols), highlightAtomColors=cols,
                                       highlightBonds=[b.GetIdx() for b in m.GetBonds() if b.GetBeginAtomIdx() in match_atoms and b.GetEndAtomIdx() in match_atoms] if match_atoms else [],
                                       highlightBondColors={b.GetIdx(): _HL for b in m.GetBonds() if b.GetBeginAtomIdx() in match_atoms and b.GetEndAtomIdx() in match_atoms})
    d.FinishDrawing(); return "data:image/png;base64," + base64.b64encode(d.GetDrawingText()).decode()


@lru_cache(maxsize=4096)
def warhead_png(ligand_sdf: str, match_atoms: tuple, reactive_atoms: tuple, w: int = 130, h: int = 90) -> str | None:
    """The warhead substructure on its own (matched atoms plus one neighbouring atom as context shown as R), reactive atom in red."""
    if not match_atoms: return None
    m = _ligand_mol(ligand_sdf); keep = set(match_atoms)
    ctx = {}
    for a in list(keep):
        for n in m.GetAtomWithIdx(a).GetNeighbors():
            if n.GetIdx() not in keep: ctx.setdefault(a, []).append(n.GetIdx())
    # fragment: matched atoms only; dangling bonds are drawn as attachment points by adding dummy atoms
    frag = Chem.RWMol()
    idx = {}
    for a in sorted(keep):
        at = m.GetAtomWithIdx(a); na = Chem.Atom(at.GetAtomicNum()); na.SetFormalCharge(at.GetFormalCharge()); na.SetIsAromatic(at.GetIsAromatic()); idx[a] = frag.AddAtom(na)
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in idx and j in idx: frag.AddBond(idx[i], idx[j], b.GetBondType())
    for a, ns in ctx.items():
        r = frag.AddAtom(Chem.Atom(0)); frag.AddBond(idx[a], r, Chem.BondType.SINGLE)
    fm = frag.GetMol()
    try: Chem.SanitizeMol(fm, catchErrors=True)
    except Exception: pass
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(fm)
    for a in fm.GetAtoms():
        if a.GetAtomicNum() == 0: a.SetProp("atomLabel", "R")
    d = rdMolDraw2D.MolDraw2DCairo(w, h); d.drawOptions().clearBackground = True
    rx = {idx[i]: _RX for i in reactive_atoms if i in idx}
    rdMolDraw2D.PrepareAndDrawMolecule(d, fm, highlightAtoms=list(rx), highlightAtomColors=rx)
    d.FinishDrawing(); return "data:image/png;base64," + base64.b64encode(d.GetDrawingText()).decode()
