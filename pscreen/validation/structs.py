"""RCSB structure utilities for validation: download (cached), split receptor / ligand, CCD SMILES, symmetry-aware RMSD."""
from __future__ import annotations
import json, urllib.request
from pathlib import Path
import numpy as np
from rdkit import Chem

CACHE = Path(__file__).resolve().parents[2] / "data" / "cache"


def _get(url: str, dest: Path) -> Path:
    if not dest.exists() or dest.stat().st_size == 0:
        dest.parent.mkdir(parents=True, exist_ok=True)
        req = urllib.request.Request(url, headers={"User-Agent": "pscreen-validation/0.1"})
        with urllib.request.urlopen(req, timeout=60) as r:
            dest.write_bytes(r.read())
    return dest


def fetch_pdb(pdb_id: str) -> Path:
    pdb_id = pdb_id.upper()
    return _get(f"https://files.rcsb.org/download/{pdb_id}.pdb", CACHE / "pdb" / f"{pdb_id}.pdb")


def ccd_smiles(comp_id: str) -> str:
    p = _get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id.upper()}", CACHE / "ccd" / f"{comp_id.upper()}.json")
    d = json.loads(p.read_text())
    desc = d.get("rcsb_chem_comp_descriptor", {})
    smi = desc.get("SMILES_stereo") or desc.get("SMILES")
    if not smi:
        raise ValueError(f"No SMILES available for CCD component {comp_id}")
    return smi


def split_complex(pdb_path: Path, receptor_chains: list[str] | None, lig_resname: str, lig_chain: str | None = None,
                  lig_resnum: int | None = None, outdir: Path | None = None) -> dict:
    """Write receptor (protein ATOM records of the chosen chains, no waters/ligands) and the ligand instance as PDB.
    Returns paths. Ligand instance = first matching HETATM group (resname[, chain][, resseq])."""
    outdir = Path(outdir or pdb_path.parent); outdir.mkdir(parents=True, exist_ok=True)
    stem = pdb_path.stem
    prot, lig = [], []
    for l in open(pdb_path):
        if l.startswith("ATOM") and (receptor_chains is None or l[21] in receptor_chains):
            if l[16] in " A": prot.append(l)
        elif l.startswith("HETATM") and l[17:20].strip() == lig_resname and (lig_chain is None or l[21] == lig_chain) \
                and (lig_resnum is None or int(l[22:26]) == lig_resnum):
            lig.append(l)
    if not lig:
        raise ValueError(f"{stem}: ligand {lig_resname} {lig_chain} {lig_resnum} not found")
    key = (lig[0][21], lig[0][22:27])   # keep a single instance
    lig = [l for l in lig if (l[21], l[22:27]) == key and l[16] in " A"]
    rec_p, lig_p = outdir / f"{stem}_rec.pdb", outdir / f"{stem}_{lig_resname}_ref.pdb"
    rec_p.write_text("".join(prot) + "END\n"); lig_p.write_text("".join(lig) + "END\n")
    return {"receptor": rec_p, "ligand_ref": lig_p, "ligand_chain": key[0], "ligand_resnum": int(key[1][:4])}


def near_chains(pdb_path: Path, lig_resname: str, lig_chain: str | None, lig_resnum: int | None, cutoff: float = 5.0) -> list[str]:
    """Protein chains with any heavy atom within `cutoff` A of the ligand instance (pockets can span subunits)."""
    lig, prot = [], []
    for l in open(pdb_path):
        if l.startswith("HETATM") and l[17:20].strip() == lig_resname and (lig_chain is None or l[21] == lig_chain) \
                and (lig_resnum is None or int(l[22:26]) == lig_resnum) and l[76:78].strip() != "H":
            lig.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif l.startswith("ATOM") and l[76:78].strip() != "H":
            prot.append((l[21], float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if not lig: raise ValueError(f"ligand {lig_resname} not found in {pdb_path.stem}")
    L = np.array(lig); P = np.array([p[1:] for p in prot]); ch = np.array([p[0] for p in prot])
    near = set()
    for i in range(0, len(P), 2000):
        d = np.linalg.norm(P[i:i+2000, None, :] - L[None, :, :], axis=2).min(1)
        near |= set(ch[i:i+2000][d <= cutoff])
    return sorted(near)


def ref_ligand_mol(ref_pdb: Path) -> Chem.Mol:
    """Crystal ligand with proximity-derived bonds (bond orders unreliable; only used as a heavy-atom graph + coordinates)."""
    m = Chem.MolFromPDBFile(str(ref_pdb), sanitize=False, removeHs=True, proximityBonding=True)
    if m is None: raise ValueError(f"cannot read {ref_pdb}")
    return m


def _skeleton(mol: Chem.Mol):
    rw = Chem.RWMol(); idx = {}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1: idx[a.GetIdx()] = rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in idx and j in idx: rw.AddBond(idx[i], idx[j], Chem.BondType.SINGLE)
    return rw.GetMol(), list(idx)


def match_graphs(ref: Chem.Mol, pose: Chem.Mol, max_matches: int = 20000) -> list[list[tuple[int, int]]]:
    """All heavy-atom element/connectivity isomorphisms ref->pose, as lists of (ref_idx, pose_idx). [] if none."""
    rs, r_orig = _skeleton(ref); ps, p_orig = _skeleton(pose)
    if rs.GetNumAtoms() != ps.GetNumAtoms(): return []
    ms = ps.GetSubstructMatches(rs, uniquify=False, maxMatches=max_matches)
    return [[(r_orig[qi], p_orig[ti]) for qi, ti in enumerate(m)] for m in ms]


def symmetric_rmsd(ref: Chem.Mol, pose: Chem.Mol) -> float | None:
    """Heavy-atom RMSD with no re-alignment (docking frame), minimised over graph automorphisms. None if graphs do not match."""
    ms = match_graphs(ref, pose)
    if not ms: return None
    rx = ref.GetConformer().GetPositions(); px = pose.GetConformer().GetPositions()
    best = min(np.sqrt(np.mean([np.sum((rx[a] - px[b]) ** 2) for a, b in m])) for m in ms)
    return float(best)
