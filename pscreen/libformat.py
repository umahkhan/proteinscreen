"""The standard compound-library format, a template for it, and a validator that explains every problem in plain words.

Standard:  a CSV with a header row and exactly these two columns
    id,smiles
    ASP-001,CC(=O)Oc1ccccc1C(=O)O
- id     unique compound name (letters, numbers, - _ . ; no spaces needed). Blank ids are numbered record_1, record_2, ...
- smiles one valid SMILES string per row (no salts needed; the largest fragment is kept).
Other inputs are still read (CSV/TSV with a smiles/canonical_smiles column, Excel, .smi, .sdf) and shown to the user converted to this form.
"""
from __future__ import annotations
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd
from .prep import read_library

TEMPLATE_CSV = ("id,smiles\n"
                "aspirin,CC(=O)Oc1ccccc1C(=O)O\n"
                "ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O\n"
                "acrylamide_example,C=CC(=O)Nc1ccccc1\n")
FORMAT_HELP = ("**Standard format:** a CSV with a header row and two columns, `id` and `smiles`, one compound per row "
               "(e.g. `aspirin,CC(=O)Oc1ccccc1C(=O)O`). Excel, TSV, `.smi` and SDF files are also accepted and converted to this form.")
_ID_OK = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-+()\[\]:/ ]{0,63}$")


@dataclass
class LibraryReport:
    table: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=["id", "smiles", "status"]))   # every input row with status 'ok' or the reason
    errors: list[str] = field(default_factory=list)       # file-level problems: nothing can be run
    warnings: list[str] = field(default_factory=list)     # row-level problems: those rows are reported and skipped, never silently dropped

    @property
    def n_total(self) -> int: return len(self.table)
    @property
    def n_ok(self) -> int: return int((self.table["status"] == "ok").sum()) if len(self.table) else 0
    @property
    def ok(self) -> bool: return not self.errors and self.n_ok > 0

    def standard_csv(self) -> str:
        """Valid rows only, in the standard `id,smiles` form."""
        t = self.table[self.table["status"] == "ok"][["id", "smiles"]]
        buf = io.StringIO(); t.to_csv(buf, index=False); return buf.getvalue()


def validate_library(path: str | Path) -> LibraryReport:
    path = Path(path); rep = LibraryReport()
    try:
        if path.suffix.lower() in (".xlsx", ".xls"):
            csv = path.with_suffix(".converted.csv"); pd.read_excel(path).to_csv(csv, index=False); path = csv
        recs = read_library(path)
    except Exception as e:
        rep.errors.append(_friendly(str(e), path)); return rep
    if not recs: rep.errors.append("The file contains no compounds."); return rep
    seen: dict[str, int] = {}; rows = []
    for r in recs:
        rows.append([str(r.original_id).strip(), r.input_smiles if r.input_smiles is not None else Chem_smiles(r.mol), r.error])
    for i, (cid, smi, err) in enumerate(rows):
        status = "ok"
        if err: status = err
        elif not _ID_OK.match(cid): status = f"id '{cid[:30]}' has unsupported characters or is too long"
        elif cid in seen: status = f"duplicate id (first seen in row {seen[cid] + 1})"
        if status == "ok": seen[cid] = i
        rows[i] = [cid, smi or "", status]
    rep.table = pd.DataFrame(rows, columns=["id", "smiles", "status"])
    bad = rep.table[rep.table["status"] != "ok"]
    if len(bad): rep.warnings.append(f"{len(bad)} of {len(rep.table)} rows will be skipped (listed below).")
    if rep.n_ok == 0: rep.errors.append("No valid compounds were found. " + FORMAT_HELP)
    return rep


def Chem_smiles(mol) -> str:
    from rdkit import Chem
    try: return Chem.MolToSmiles(mol) if mol is not None else ""
    except Exception: return ""


def _friendly(msg: str, path: Path) -> str:
    if "No SMILES column" in msg: return "Could not find a SMILES column. Name the columns `id` and `smiles` (see the template). " + msg.split(";")[-1].strip()
    if "Unsupported library format" in msg: return f"'{path.suffix}' files are not supported. Use CSV (id,smiles), Excel, TSV, .smi or SDF."
    return f"Could not read the compound file: {msg}"
