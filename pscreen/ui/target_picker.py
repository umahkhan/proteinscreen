"""Target-residue chooser with validation against the structure. Used on New screen and in the Review sidebar."""
from __future__ import annotations
import streamlit as st
from ..residues import REACTIVE_TYPES, near_pocket, residue_table, validate_target


def pick(receptor_pdb: str, key: str, default: str = "", pocket_center=None, container=None) -> tuple[str, dict | None]:
    """Returns (validated spec or '', validation dict or None). Never returns an unvalidated residue."""
    c = container or st
    try: table = residue_table(receptor_pdb)
    except Exception as e: c.error(f"Could not read the structure: {e}"); return "", None
    mode = c.radio("Target residue", ["None", "Pick from structure", "Type it"], horizontal=True, key=f"{key}_mode", index=(1 if default else 0),
                   help="The residue whose nucleophile should reach the ligand's warhead. Enables warhead-to-residue geometry and the covalent lane.")
    if mode == "None": return "", None
    spec = ""
    if mode == "Pick from structure":
        rt = c.selectbox("Residue type", REACTIVE_TYPES, key=f"{key}_type", index=REACTIVE_TYPES.index(default.split(":")[1][:3]) if default and default.split(":")[1][:3] in REACTIVE_TYPES else 0)
        sub = table[(table.resname == rt) & table.nuc_xyz.notna()]
        if sub.empty: c.warning(f"No {rt} residues with a complete side chain in this structure."); return "", None
        if pocket_center is not None:
            nr = near_pocket(table, pocket_center, 1e6, [rt]); sub = sub.merge(nr[["chain", "resnum", "dist_to_pocket"]], on=["chain", "resnum"]).sort_values("dist_to_pocket")
        opts = {f"{r.chain}:{r.resname}{r.resnum}" + (f"  ·  {r.dist_to_pocket} Å from pocket" if pocket_center is not None else ""): f"{r.chain}:{r.resname}{r.resnum}" for r in sub.itertuples()}
        labels = list(opts)
        idx = next((i for i, l in enumerate(labels) if opts[l] == default), 0)
        spec = opts[c.selectbox(f"{rt} residue ({len(labels)} in structure" + (", nearest to the pocket first)" if pocket_center is not None else ")"), labels, index=idx, key=f"{key}_res")]
    else:
        spec = c.text_input("Residue", default, placeholder="A:CYS145", key=f"{key}_txt", help="CHAIN:RESNAMEnumber. Also accepted: CYS145 (single-chain), A:145.")
        if not spec: return "", None
    v = validate_target(receptor_pdb, spec, pocket_center, table=table)
    if not v["ok"]:
        for e in v["errors"]: c.error(e)
        return "", v
    c.success(f"{v['spec']} found · nucleophile {v['nucleophile_atom']}" + (f" · {v['distance_to_pocket']} Å from pocket centre" if v["distance_to_pocket"] is not None else ""))
    for w in v["warnings"]: c.warning(w)
    return v["spec"], v
