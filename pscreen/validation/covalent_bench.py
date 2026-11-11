"""Covalent pose-recovery benchmark on the GalaxyCDock covalent set (seoklab/cdock_dataset, cov_set.json).

Per entry we know the bonded protein residue/atom and the ligand atom name. The ligand input is the CCD component SMILES, whose
atoms are exactly the deposited atoms (i.e. already in bound-state chemistry), so this isolates GNINA's covalent docking from our
warhead-transform library (which is unit-tested separately). Each ligand is docked twice into the same receptor: covalently and
noncovalently. Success = symmetric heavy-atom RMSD <= 2 A vs the crystal pose.
This measures POSE recovery only. It says nothing about reaction kinetics, potency or selectivity.
Split: at most one complex per UniProt target (target-level de-duplication). Leakage caveat: GNINA's CNN scorers were trained on
PDBbind/CrossDocked-style data (mostly noncovalent) whose overlap with this set is unknown.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
from ..docking import DockConfig, run_screen, MODE_COVALENT, MODE_NONCOVALENT
from ..prep import prepare_library
from ..pockets import from_known_ligand
from ..covalent import attach_smarts, bound_forms
from ..warheads import find_warheads
from .structs import fetch_pdb, ccd_smiles, split_complex, ref_ligand_mol, match_graphs, _skeleton, near_chains

COV_SET = Path(__file__).resolve().parents[2] / "data" / "cdock" / "test_data" / "cov_set.json"


def select_entries(n: int = 30, seed: int = 0, residues=("CYS", "SER", "LYS", "TYR", "HIS"), max_res: float = 2.5,
                   heavy=(10, 45), max_tor: int = 10, per_residue_cap: dict | None = None) -> list[dict]:
    d = [e for e in json.load(open(COV_SET)) if e["exp_type"] == "X-RAY DIFFRACTION" and e["res_name"] in residues
         and isinstance(e["resolution"], (int, float)) and e["resolution"] <= max_res
         and heavy[0] <= e["lig_atom_num"] <= heavy[1] and e["lig_n_tor"] <= max_tor and str(e["lig_num_auth"]).lstrip("-").isdigit()]
    best: dict[str, dict] = {}
    for e in sorted(d, key=lambda e: (e["resolution"], e["pdb"])):         # one per UniProt target (best resolution)
        best.setdefault(e["uniprot_id"], e)
    pool = list(best.values())
    rng = np.random.default_rng(seed); rng.shuffle(pool)
    cap = per_residue_cap or {"CYS": int(n * 0.6), "SER": int(n * 0.2), "LYS": int(n * 0.1), "TYR": int(n * 0.05), "HIS": int(n * 0.05)}
    out, cnt = [], {}
    for e in pool:
        if cnt.get(e["res_name"], 0) < cap.get(e["res_name"], 0):
            out.append(e); cnt[e["res_name"]] = cnt.get(e["res_name"], 0) + 1
        if len(out) >= n: break
    return out


def build(e: dict, workdir: Path, form: str = "ccd") -> dict:
    """form='ccd': dock the CCD SMILES as deposited. form='pscreen': if a recognised warhead covers the bonded atom and the adduct keeps
    the deposited heavy-atom graph, dock our adduct form instead (else falls back to 'ccd'; recorded in `ligand_form`)."""
    key = f"{e['pdb'].upper()}_{e['lig_name']}"
    pdb = fetch_pdb(e["pdb"])
    chains = sorted(set(near_chains(pdb, e["lig_name"], e["lig_chain_auth"], int(e["lig_num_auth"]))) | {e["res_chain_auth"]})
    parts = split_complex(pdb, chains, e["lig_name"], e["lig_chain_auth"], int(e["lig_num_auth"]), workdir / key)
    ref = ref_ligand_mol(parts["ligand_ref"])
    smi = ccd_smiles(e["lig_name"]); ccd = Chem.MolFromSmiles(smi)
    if ccd is None: raise ValueError(f"CCD SMILES unparsable: {smi}")
    ridx = next((a.GetIdx() for a in ref.GetAtoms() if a.GetPDBResidueInfo().GetName().strip() == e["lig_atom"]), None)
    if ridx is None: raise ValueError(f"ligand atom {e['lig_atom']} not found in crystal ligand")
    # Map deposited atoms onto the CCD (free) ligand. The deposited ligand can lack leaving atoms (F-, H2O ...) present in the CCD.
    ms = match_graphs(ref, ccd)
    if ms:
        m0 = dict(ms[0])
    else:
        qs, q_orig = _skeleton(ref); ts, t_orig = _skeleton(ccd)
        sub = ts.GetSubstructMatch(qs) if qs.GetNumAtoms() < ts.GetNumAtoms() else ()
        if not sub: raise ValueError("deposited ligand atoms are not a subgraph of the CCD ligand")
        m0 = {q_orig[qi]: t_orig[ti] for qi, ti in enumerate(sub)}
    keep = set(m0.values()); lost = [a.GetSymbol() for a in ccd.GetAtoms() if a.GetIdx() not in keep]
    sidx0 = m0[ridx]
    # Form A: CCD graph restricted to the deposited atoms (CCD bond orders)
    rw = Chem.RWMol(ccd); Chem.Kekulize(rw, clearAromaticFlags=True)
    for i in sorted((a.GetIdx() for a in ccd.GetAtoms() if a.GetIdx() not in keep), reverse=True): rw.RemoveAtom(i)
    order = sorted(keep); remap = {old: new for new, old in enumerate(order)}
    for a_ in rw.GetAtoms(): a_.SetNoImplicit(False); a_.SetNumExplicitHs(0)
    try: Chem.SanitizeMol(rw)
    except Exception as ex: raise ValueError(f"cannot sanitize deposited-graph ligand: {ex}")
    smol, sidx, ligand_form = rw.GetMol(), remap[sidx0], "deposited_graph" + (f" (removed {'+'.join(lost)})" if lost else "")
    # Form B: pscreen adduct transform, only if it reproduces the deposited heavy-atom graph
    if form == "pscreen":
        for h in find_warheads(ccd):
            if sidx0 not in h["reactive_atoms"]: continue
            for bf in bound_forms(ccd, h):
                if bf.attach_idx == bf.source_to_product.get(sidx0) and bf.product.GetNumHeavyAtoms() == ref.GetNumHeavyAtoms() and match_graphs(ref, bf.product):
                    smol, sidx, ligand_form = bf.product, bf.attach_idx, f"pscreen_adduct:{h['warhead']}"
                    break
            if ligand_form.startswith("pscreen_adduct"): break
    smi = Chem.MolToSmiles(smol)
    if not match_graphs(ref, smol): raise ValueError("internal: prepared ligand graph does not match deposited ligand")
    # receptor atom coordinate (for bond-length checks)
    rxyz = None
    for l in open(parts["receptor"]):
        if l.startswith("ATOM") and l[21] == e["res_chain_auth"] and int(l[22:26]) == int(e["res_num_auth"]) and l[12:16].strip() == e["res_atom"]:
            rxyz = (float(l[30:38]), float(l[38:46]), float(l[46:54])); break
    if rxyz is None: raise ValueError(f"receptor atom {e['res_chain_auth']}:{e['res_num_auth']}:{e['res_atom']} not found")
    cx = ref.GetConformer().GetPositions()[ridx]
    return {**e, "key": key, "smiles": smi, "ligand_form": ligand_form, "receptor": str(parts["receptor"]), "ligand_ref": str(parts["ligand_ref"]),
            "rec_atom": f"{e['res_chain_auth']}:{int(e['res_num_auth'])}:{e['res_atom']}", "lig_smarts": attach_smarts(smol, sidx),
            "ref_idx": ridx, "rec_xyz": rxyz, "crystal_bond_A": round(float(np.linalg.norm(cx - np.array(rxyz))), 2)}


def run(entries: list[dict], outdir: str | Path, cfg: DockConfig, progress=None, form: str = "ccd",
        modes: tuple = (MODE_COVALENT, MODE_NONCOVALENT), only_adducts: bool = False) -> pd.DataFrame:
    """form 'ccd' = deposited-graph ligand; 'pscreen' = our adduct transform where applicable. only_adducts keeps just the entries
    where the transform changed the ligand (used to compare forms on exactly those complexes)."""
    outdir = Path(outdir).resolve(); outdir.mkdir(parents=True, exist_ok=True)
    tasks, pockets, rows = [], [], []
    for e in entries:
        try: tasks.append(build(e, outdir / "inputs", form))
        except Exception as ex: rows.append({"key": f"{e['pdb']}_{e['lig_name']}", "status": "input_failed", "reason": str(ex), "residue": e["res_name"]})
    if only_adducts:
        tasks = [t for t in tasks if t["ligand_form"].startswith("pscreen_adduct")]
        if not tasks: raise ValueError("no entries where the pscreen adduct transform applies")
    for t in tasks:
        pk = from_known_ligand(f"{t['key']}_site", t["receptor"], t["ligand_ref"], padding=6.0, min_size=20.0); pk.receptor_path = t["receptor"]
        pockets.append(pk)
    pd.DataFrame([{"id": t["key"], "smiles": t["smiles"]} for t in tasks]).to_csv(outdir / "ligands.csv", index=False)
    prepared = prepare_library(outdir / "ligands.csv", outdir)
    pairs = {(t["key"], f"{t['key']}_site") for t in tasks}
    specs = {(t["key"], f"{t['key']}_site"): {"rec_atom": t["rec_atom"], "lig_smarts": t["lig_smarts"], "bond_order": 1} for t in tasks}
    store = None
    for mode in modes:
        store = run_screen(prepared, pockets[0].receptor_path, pockets, cfg, outdir, progress=progress, pairs=pairs, mode=mode, covalent_specs=specs)
    return score(outdir, tasks, prepared, store, rows, modes)


def score(outdir: Path, tasks, prepared, store, rows, modes=(MODE_COVALENT, MODE_NONCOVALENT)) -> pd.DataFrame:
    poses, jobs = store.poses(), store.jobs()
    for t in tasks:
        base = {"key": t["key"], "pdb": t["pdb"], "residue": t["res_name"], "rec_atom": t["rec_atom"], "crystal_bond_A": t["crystal_bond_A"],
                "heavy_atoms": t["lig_atom_num"], "ligand_form": t["ligand_form"], "resolution": t["resolution"], "uniprot": t["uniprot_id"]}
        prow = prepared[prepared.ligand_key == t["key"]]
        if prow.empty or prow.iloc[0].status != "ok":
            rows.append({**base, "status": "prep_failed", "reason": None if prow.empty else prow.iloc[0].reason}); continue
        ref = ref_ligand_mol(t["ligand_ref"]); rec = np.array(t["rec_xyz"])
        for mode in modes:
            g = poses[(poses.ligand_key == t["key"]) & (poses["mode"] == mode)]
            r = {**base, "mode": mode, "status": "ok", "reason": None}
            if g.empty:
                j = jobs[(jobs.ligand_key == t["key"]) & (jobs["mode"] == mode)]
                r.update(status="docking_failed", reason=j.iloc[0].error if len(j) else "not run"); rows.append(r); continue
            sup = [m for m in Chem.SDMolSupplier(g.iloc[0].pose_sdf, removeHs=False, sanitize=False)]
            rm, bond = [], []
            for i in g.pose_idx:
                pm = Chem.RemoveHs(sup[int(i)], sanitize=False)
                ms = match_graphs(ref, pm)
                if not ms: rm.append(np.nan); bond.append(np.nan); continue
                px, rx = pm.GetConformer().GetPositions(), ref.GetConformer().GetPositions()
                best = min(ms, key=lambda m: np.mean([np.sum((rx[a] - px[b]) ** 2) for a, b in m]))
                rm.append(float(np.sqrt(np.mean([np.sum((rx[a] - px[b]) ** 2) for a, b in best]))))
                bond.append(float(np.linalg.norm(px[dict(best)[t["ref_idx"]]] - rec)))
            g = g.assign(rmsd=rm, bond=bond)
            if g.rmsd.isna().all(): r.update(status="rmsd_unmapped", reason="docked ligand graph does not match crystal ligand"); rows.append(r); continue
            for name, col, asc in [("vina_like", "minimized_affinity", True), ("cnn_pose", "cnn_score", False), ("cnn_vs", "cnn_vs", False)]:
                top = g.sort_values(col, ascending=asc).iloc[0]; r[f"rmsd_top1_{name}"] = round(float(top.rmsd), 2)
            top = g.sort_values("cnn_score", ascending=False).iloc[0]; r["attach_to_rec_atom_A_top1_cnn_pose"] = round(float(top.bond), 2)
            r["n_poses_bond_satisfied"] = int((g.bond <= 2.3).sum())
            gb = g[g.bond <= 2.3]       # principled covalent selection: only poses in which the modelled bond is formed
            for name, col, asc in [("vina_like", "minimized_affinity", True), ("cnn_pose", "cnn_score", False)]:
                r[f"rmsd_top1_{name}_bondfiltered"] = round(float(gb.sort_values(col, ascending=asc).iloc[0].rmsd), 2) if len(gb) else np.nan
            r["rmsd_best_any_pose_oracle"] = round(float(g.rmsd.min()), 2); r["n_poses"] = len(g)
            rows.append(r)
    res = pd.DataFrame(rows); res.to_csv(outdir / "covalent_results.csv", index=False)
    return res


def summarize(res: pd.DataFrame, cutoff: float = 2.0) -> pd.DataFrame:
    ok = res[res.status == "ok"]
    cols = [c for c in ok.columns if c.startswith("rmsd_")]
    out = ok.groupby("mode")[cols].agg(lambda s: round(float((s <= cutoff).mean()), 3)).add_prefix("success_")
    out["n_scored"] = ok.groupby("mode").size()
    return out
