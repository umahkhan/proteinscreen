"""Target-residue discovery and validation against the uploaded structure.

A target residue spec is 'CHAIN:RESNAMEnumber' (e.g. A:CYS145). Validation checks that the chain and residue exist, that the residue name
matches (a user asking for 'cysteine 145' on a structure where 145 is a serine gets a clear error), that the nucleophilic atom exists in the
file, and, when a pocket centre is known, how far the residue is from the pocket.
"""
from __future__ import annotations
import re
import numpy as np, pandas as pd
from .geometry import NUCLEOPHILE_ATOMS, load_structure

THREE = {"C": "CYS", "S": "SER", "K": "LYS", "Y": "TYR", "H": "HIS", "T": "THR", "D": "ASP", "E": "GLU"}
REACTIVE_TYPES = ["CYS", "SER", "LYS", "TYR", "HIS", "THR"]


def residue_table(pdb_path: str) -> pd.DataFrame:
    """One row per standard residue: chain, resnum, resname, nucleophile atom (first available), its xyz, CA xyz."""
    rows = []
    for model in load_structure(pdb_path):
        for ch in model:
            for r in ch:
                if r.id[0].strip() != "": continue
                names = NUCLEOPHILE_ATOMS.get(r.get_resname(), [])
                nuc = next((n for n in names if n in r), None)
                rows.append({"chain": ch.id, "resnum": r.id[1], "icode": r.id[2].strip(), "resname": r.get_resname(), "nucleophile_atom": nuc,
                             "nuc_xyz": tuple(r[nuc].coord.astype(float)) if nuc else None, "ca_xyz": tuple(r["CA"].coord.astype(float)) if "CA" in r else None})
        break
    return pd.DataFrame(rows)


def parse_spec(spec: str) -> tuple[str | None, int | None, str | None, str | None]:
    """Accepts 'A:CYS145', 'A:145', 'CYS145', 'Cys 145', 'A145'. Returns (chain, num, resname, error)."""
    s = (spec or "").strip().upper().replace(" ", "")
    if not s: return None, None, None, "No residue given."
    m = re.fullmatch(r"(?:([A-Z0-9]):)?([A-Z]{3})?(-?\d+)", s) or re.fullmatch(r"([A-Z0-9])(-?\d+)", s)
    if not m: return None, None, None, f"Could not read '{spec}'. Use CHAIN:RESNAMEnumber, e.g. A:CYS145."
    if m.re.pattern.startswith("(?:"): chain, name, num = m.group(1), m.group(2), int(m.group(3))
    else: chain, name, num = m.group(1), None, int(m.group(2))
    return chain, num, name, None


def validate_target(pdb_path: str, spec: str, pocket_center=None, max_pocket_dist: float = 15.0, table: pd.DataFrame | None = None) -> dict:
    """Returns {ok, errors[], warnings[], spec (normalised), resname, nucleophile_atom, distance_to_pocket}."""
    out = {"ok": False, "errors": [], "warnings": [], "spec": None, "resname": None, "nucleophile_atom": None, "distance_to_pocket": None}
    chain, num, name, err = parse_spec(spec)
    if err: out["errors"].append(err); return out
    t = table if table is not None else residue_table(pdb_path)
    if t.empty: out["errors"].append("The structure contains no standard protein residues."); return out
    chains = sorted(t.chain.unique())
    if chain is None:
        if len(chains) > 1: out["errors"].append(f"Chain not given and the structure has chains {chains}. Write it as CHAIN:RESNAMEnumber."); return out
        chain = chains[0]
    if chain not in chains: out["errors"].append(f"Chain '{chain}' is not in the structure (chains: {', '.join(chains)})."); return out
    sub = t[(t.chain == chain) & (t.resnum == num)]
    if sub.empty:
        same = t[(t.chain == chain) & (t.resname == name)] if name else t.iloc[0:0]
        hint = f" {name} residues on chain {chain}: {', '.join(f'{r.resname}{r.resnum}' for r in same.itertuples())}." if len(same) else ""
        lo, hi = t[t.chain == chain].resnum.min(), t[t.chain == chain].resnum.max()
        out["errors"].append(f"Residue {num} does not exist on chain {chain} (numbering {lo}–{hi}).{hint}"); return out
    r = sub.iloc[0]
    if name and r.resname != name:
        near = t[(t.chain == chain) & (t.resname == name)]
        hint = f" {name} residues on chain {chain}: {', '.join(str(x) for x in near.resnum)}." if len(near) else f" Chain {chain} has no {name}."
        out["errors"].append(f"Residue {num} on chain {chain} is {r.resname}, not {name}.{hint}"); return out
    out.update(spec=f"{chain}:{r.resname}{num}", resname=r.resname, nucleophile_atom=r.nucleophile_atom)
    if r.resname not in NUCLEOPHILE_ATOMS: out["errors"].append(f"{r.resname} has no nucleophilic side-chain atom supported for covalent analysis (supported: {', '.join(NUCLEOPHILE_ATOMS)})."); return out
    if r.nucleophile_atom is None: out["errors"].append(f"{r.resname}{num} is missing its nucleophilic side-chain atom ({'/'.join(NUCLEOPHILE_ATOMS[r.resname])}) in the file (truncated side chain?)."); return out
    if r.icode: out["warnings"].append(f"Insertion code '{r.icode}' present; residue identified by number only.")
    if pocket_center is not None:
        d = float(np.linalg.norm(np.array(r.nuc_xyz) - np.array(pocket_center))); out["distance_to_pocket"] = round(d, 1)
        if d > max_pocket_dist: out["warnings"].append(f"{out['spec']} is {d:.1f} Å from the pocket centre; its warhead geometry will be far from the docked ligands.")
    out["ok"] = True
    return out


def near_pocket(table: pd.DataFrame, center, radius: float = 15.0, names=REACTIVE_TYPES) -> pd.DataFrame:
    t = table[table.resname.isin(names) & table.nuc_xyz.notna()].copy()
    t["dist_to_pocket"] = [round(float(np.linalg.norm(np.array(x) - np.array(center))), 1) for x in t.nuc_xyz]
    return t[t.dist_to_pocket <= radius].sort_values("dist_to_pocket")
