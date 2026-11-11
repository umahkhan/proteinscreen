"""Pockets-only helpers: fetch a structure, summarise pockets as a table, build the download bundle. No docking involved."""
from __future__ import annotations
import io, json, re, urllib.request, zipfile
from pathlib import Path
import pandas as pd
from . import pymol_export

LIMITS = {"max_mb": 15, "max_residues": 2500, "lacuna_max_residues": 1500}
MODEL_OF = {"p2rank": "P2Rank", "pocketminer": "PocketMiner", "lacuna": "Lacuna"}


def fetch_pdb(pdb_id: str, dest: Path) -> Path:
    """Download an entry from the RCSB PDB (PDB format, falling back to mmCIF)."""
    pid = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pid): raise ValueError("A PDB ID is 4 characters, e.g. 3MXF.")
    for ext in ("pdb", "cif"):
        try:
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.{ext}", timeout=30) as r: data = r.read()
            if len(data) > LIMITS["max_mb"] * 2**20: raise ValueError("Entry is larger than the size limit.")
            out = dest / f"{pid}.{ext}"; out.write_bytes(data); return out
        except urllib.error.HTTPError as e:
            if e.code != 404: raise
    raise ValueError(f"No entry {pid} found at the RCSB PDB.")


def model_of(p) -> str:
    return MODEL_OF.get(p.pocket_id.split("_")[0], p.pocket_id.split("_")[0])


def pocket_dict(p) -> dict:
    return {"id": p.pocket_id, "residues": p.residues, "center": list(p.center), "size": list(p.size), "source": p.source,
            "residue_scores": p.metrics.get("residue_scores"), "residue_score_label": p.metrics.get("residue_score_label", "")}


def pockets_table(pockets) -> pd.DataFrame:
    rows = []
    for p in pockets:
        m = p.metrics or {}
        score = m.get("probability") if m.get("probability") is not None else m.get("persistence")
        rows.append({"pocket": p.pocket_id, "tool": model_of(p), "type": "cryptic (predicted)" if p.source == "predicted_cryptic" else "conventional (predicted)",
                     "tool score": m.get("score"), "probability / persistence": score, "volume Å³": m.get("volume_open_A3"), "residues": len(p.residues),
                     "center x,y,z": ", ".join(f"{c:.1f}" for c in p.center), "lining residues": " ".join(p.residues)})
    return pd.DataFrame(rows)


README = """pscreen pockets export

pockets.csv     one row per predicted pocket (tool, type, score, lining residues, centre)
pockets.json    the same with box sizes and tool metrics
receptor.pdb    the structure the tools were run on (waters, ions and ligands removed)
pymol_pockets/  PyMOL scene (unzip, then in PyMOL:  @pscreen_scene.pml)

Pockets are computational predictions. They are not experimentally confirmed and do not show that a compound binds, or that a pocket is druggable.
Cryptic-pocket predictions on closed (apo) structures are hypotheses. Tools: P2Rank (MIT), PocketMiner (MIT), Lacuna (MIT); please cite them.
"""


def bundle(receptor_path: str, pockets, status: dict, panel: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("pscreen_pockets/README.txt", README)
        z.writestr("pscreen_pockets/pockets.csv", pockets_table(pockets).to_csv(index=False))
        z.writestr("pscreen_pockets/pockets.json", json.dumps({"tool_status": status, "pockets": [
            {**pocket_dict(p), "metrics": {k: v for k, v in (p.metrics or {}).items() if k != "residue_scores"}} for p in pockets]}, indent=1, default=str))
        z.writestr("pscreen_pockets/receptor.pdb", Path(receptor_path).read_text())
        for name, text in pymol_export.build(panel).items(): z.writestr(f"pscreen_pockets/pymol_pockets/{name}", text)
    return buf.getvalue()
