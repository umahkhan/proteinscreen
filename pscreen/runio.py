"""Load a finished or in-progress run directory into the objects the UI/analysis need."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd
from .docking import Store, Pocket
from .results import build_results


def _csv(path: Path):
    """Optional run files can be missing or still empty while the pipeline is writing them."""
    try: return pd.read_csv(path) if path.exists() and path.stat().st_size > 0 else None
    except (pd.errors.EmptyDataError, pd.errors.ParserError): return None


def list_runs(root="runs") -> list[Path]:
    return sorted((p.parent for p in Path(root).glob("*/run_config.json")), key=lambda p: p.stat().st_mtime, reverse=True)


def load_pockets(d: list[dict]) -> list[Pocket]:
    return [Pocket(**{**x, "center": tuple(x["center"]), "size": tuple(x["size"])}) for x in d]


def load_run(rundir: Path, target: str | None = None) -> dict:
    """target: optional override of the target residue used for noncovalent geometry columns (covalent lane stays as docked)."""
    store = Store(rundir / "results.sqlite")
    meta = {k: json.loads(v) for k, v in store.db.execute("SELECT key,value FROM meta")}
    pockets = load_pockets(meta["pockets"])
    extra = Path(rundir) / "pockets_display.json"      # overlay-only predictions (no docking at these sites)
    if extra.exists(): pockets += load_pockets(json.loads(extra.read_text()))
    prepared = pd.read_csv(rundir / "prepared.csv")
    cmap = _csv(rundir / "covalent_map.csv")
    props = _csv(rundir / "props" / "admet_out.csv")
    if props is not None:
        eps = [c for c in props.columns if c != "ligand_key"]; props = props.rename(columns={e: f"chemprop__{e}" for e in eps})
    cfg = json.loads((rundir / "run_config.json").read_text())
    exp = pd.read_csv(cfg["experimental_csv"]) if cfg.get("experimental_csv") and Path(cfg["experimental_csv"]).exists() else None
    receptor = str(rundir / "inputs" / "receptor.pdb") if (rundir / "inputs" / "receptor.pdb").exists() else meta["receptor"]
    tgt = target if target is not None else (meta.get("target_residue") or None)
    res = build_results(prepared, store, pockets, receptor, tgt or None, exp, cmap, props)
    manifest = json.loads((rundir / "run_manifest.json").read_text()) if (rundir / "run_manifest.json").exists() else {}
    admet_meta = json.loads((rundir / "props" / "admet_models.json").read_text()) if (rundir / "props" / "admet_models.json").exists() else None
    return {"store": store, "meta": meta, "pockets": pockets, "prepared": prepared, "covalent_map": cmap, "results": res, "receptor": receptor,
            "target": (tgt or ""), "docked_target": meta.get("target_residue") or "", "manifest": manifest, "experimental": exp, "config": cfg, "admet_meta": admet_meta,
            "covalent_report": _csv(rundir / "covalent_report.csv")}
