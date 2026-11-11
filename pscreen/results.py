"""Assemble the per-pose results table. Evidence signals stay in separate columns; there is no composite score.

Docking modes are kept apart: 'noncovalent' (ordinary GNINA docking of the intact compound) and 'covalent' (GNINA covalent docking of the
bound-state adduct form, see covalent.py). Their scores are NOT on a common scale and are never ranked against each other.
"""
from __future__ import annotations
import json
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from .docking import Store, Pocket, MODE_COVALENT, MODE_NONCOVALENT
from .geometry import load_structure, residue_atoms, reactive_geometry

DESC = ["mol_weight", "clogp", "tpsa", "hbd", "hba", "rot_bonds", "heavy_atoms"]
BOND_MAX_A = 2.3     # covalent bond ~1.8 A (C-S); tolerance for a 'satisfied' modelled bond


def build_results(prepared: pd.DataFrame, store: Store, pockets: list[Pocket], receptor: str, target_residue: str | None = None,
                  experimental: pd.DataFrame | None = None, covalent_map: pd.DataFrame | None = None,
                  props: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per pose; compounds with no poses (invalid/duplicate/failed) get one row carrying their status.
    covalent_map: rows (cov_key, ligand_key, warhead, attach_idx, source_to_product json) linking covalent jobs to parent compounds.
    props: optional Chemprop/ADMET table (ligand_key + chemprop__* columns), joined unchanged."""
    poses, jobs = store.poses(), store.jobs()
    pk = {p.pocket_id: p for p in pockets}
    targets: dict = {}
    def target_for(pocket_id, residue=None):
        residue = residue or target_residue
        rp = pk[pocket_id].receptor_path or receptor
        if residue and (rp, residue) not in targets:
            try: targets[(rp, residue)] = residue_atoms(load_structure(rp), residue)
            except ValueError: targets[(rp, residue)] = None
        return targets.get((rp, residue))
    prep = prepared.set_index("ligand_key", drop=False)
    cmap = covalent_map.set_index("cov_key") if covalent_map is not None and len(covalent_map) else None
    rows, sdf_cache, lig_cache = [], {}, {}
    for _, p in poses.iterrows():
        is_cov = p["mode"] == MODE_COVALENT
        parent_key = cmap.loc[p.ligand_key, "ligand_key"] if is_cov and cmap is not None and p.ligand_key in cmap.index else p.ligand_key
        meta = prep.loc[parent_key]
        row = {"compound_id": meta.compound_id, "ligand_key": parent_key, "pocket_id": p.pocket_id, "pocket_source": pk[p.pocket_id].source,
               "pocket_evidence": "experimental_ligand_site" if pk[p.pocket_id].experimentally_confirmed else "computational_prediction_or_user_box",
               "pose_id": f"{p.ligand_key}|{p.pocket_id}|{p['mode']}|{p.pose_idx}", "pose_idx": p.pose_idx, "docking_mode": p["mode"],
               "gnina_minimized_affinity_kcal": p.minimized_affinity, "gnina_cnn_pose_score": p.cnn_score, "gnina_cnn_affinity_pK": p.cnn_affinity,
               "gnina_cnn_vs": p.cnn_vs, "gnina_cnn_affinity_variance": p.cnn_affinity_variance, "pose_sdf": p.pose_sdf,
               "warhead_classes": meta.get("warhead_classes") or "", "reactive_atom_indices": meta.get("reactive_atom_indices"),
               "canonical_smiles": meta.canonical_smiles, **{d: meta.get(d) for d in DESC}, "processing_status": "docked", "warnings": meta.get("warnings") or "",
               "target_residue": target_residue, "reactive_to_target_distance": None, "nucleophile_atom": None, "covalent_geometry": "not_assessed",
               "covalent_docking_status": "covalent_docking" if is_cov else "noncovalent_docking", "geometry_warhead": None}
        cres = (cmap.loc[p.ligand_key].get("residue") if is_cov and cmap is not None and p.ligand_key in cmap.index else None)
        cres = cres if isinstance(cres, str) and cres else None
        target = target_for(p.pocket_id, cres) if (target_residue or cres) else None
        if cres: row["target_residue"] = cres
        if is_cov and cmap is not None and p.ligand_key in cmap.index:
            c = cmap.loc[p.ligand_key]; row["geometry_warhead"] = c.warhead; row["warhead_classes"] = c.warhead
            row["covalent_adduct_smiles"] = c.get("adduct_smiles")
            if target is not None:
                # distance from the modelled attachment atom to the nucleophile = the covalent bond length GNINA built
                if p.pose_sdf not in sdf_cache:
                    sdf_cache[p.pose_sdf] = [m for m in Chem.SDMolSupplier(p.pose_sdf, removeHs=False, sanitize=False)]
                pose = sdf_cache[p.pose_sdf][int(p.pose_idx)]
                lig = Chem.RemoveHs(Chem.SDMolSupplier(c.adduct_sdf, removeHs=False)[0])
                g = reactive_geometry(lig, pose, [int(c.attach_idx)], target)
                label = "covalent_bond_modelled" if g.distance is not None and g.distance <= BOND_MAX_A else ("bond_not_formed" if g.distance is not None else g.label)
                row.update(reactive_to_target_distance=g.distance, nucleophile_atom=g.nucleophile_atom, covalent_geometry=label)
        else:
            hits = json.loads(meta.warheads) if isinstance(meta.get("warheads"), str) else []
            if target is not None and hits:
                if p.pose_sdf not in sdf_cache:
                    sdf_cache[p.pose_sdf] = [m for m in Chem.SDMolSupplier(p.pose_sdf, removeHs=False, sanitize=False)]
                if meta.ligand_sdf not in lig_cache:
                    lig_cache[meta.ligand_sdf] = Chem.RemoveHs(Chem.SDMolSupplier(meta.ligand_sdf, removeHs=False)[0])
                pose = sdf_cache[p.pose_sdf][int(p.pose_idx)]
                best = None
                for h in hits:
                    g = reactive_geometry(lig_cache[meta.ligand_sdf], pose, h["reactive_atoms"], target)
                    if g.distance is not None and (best is None or g.distance < best[0].distance): best = (g, h["warhead"])
                    elif best is None and g.label == "unmapped": best = (g, h["warhead"])
                if best:
                    g, wname = best
                    row.update(covalent_docking_status="noncovalent_pose_assessed_for_covalent_geometry", reactive_to_target_distance=g.distance,
                               nucleophile_atom=g.nucleophile_atom, covalent_geometry=g.label, geometry_warhead=wname)
        rows.append(row)
    res = pd.DataFrame(rows)
    # GNINA covalent mode also returns poses in which the covalent bond is NOT satisfied (the CNN rescoring does not enforce it).
    # Only poses with attachment atom <= BOND_MAX_A from the nucleophile are valid covalent poses; flag compounds that have none.
    if len(res) and (res.docking_mode == MODE_COVALENT).any():
        c = res[res.docking_mode == MODE_COVALENT]
        valid = set(c[c.covalent_geometry == "covalent_bond_modelled"].ligand_key); tried = set(c.ligand_key)
        res["covalent_note"] = res.ligand_key.map(lambda k: "covalent docking produced no pose with the bond satisfied" if k in tried and k not in valid else None)
    have = set(res["ligand_key"]) if len(res) else set()
    fail = jobs[jobs.status == "failed"]
    fail_by_parent = {}
    for _, f in fail.iterrows():
        parent = cmap.loc[f.ligand_key, "ligand_key"] if f["mode"] == MODE_COVALENT and cmap is not None and f.ligand_key in cmap.index else f.ligand_key
        fail_by_parent.setdefault(parent, []).append(f"{f['mode']} @ {f.pocket_id}: {f.error}")
    extra = []
    for _, m in prepared.iterrows():
        if m.ligand_key in have: continue
        if m.status != "ok": st, why = m.status, m.reason
        elif m.ligand_key in fail_by_parent: st, why = "docking_failed", " | ".join(fail_by_parent[m.ligand_key])
        else: st, why = "not_run", None
        extra.append({"compound_id": m.compound_id, "ligand_key": m.ligand_key, "processing_status": st, "failure_reason": why, "warnings": m.get("warnings") or "",
                      "warhead_classes": m.get("warhead_classes") or "", "canonical_smiles": m.get("canonical_smiles")})
    res = pd.concat([res, pd.DataFrame(extra)], ignore_index=True)
    # partial failures (e.g. covalent job failed but noncovalent succeeded) are surfaced as a warning on the docked rows
    if fail_by_parent and len(res):
        res["partial_failures"] = res.ligand_key.map(lambda k: " | ".join(fail_by_parent.get(k, [])) or None)
    if props is not None and len(props):
        res = res.merge(props, on="ligand_key", how="left")
    if experimental is not None:
        ex = experimental.rename(columns={c: f"exp_{c}" for c in experimental.columns if c != "compound_id"})
        res = res.merge(ex, on="compound_id", how="left")
    return res


def best_pose_per_compound(res: pd.DataFrame, sort_by: str = "gnina_cnn_vs", mode: str | None = MODE_NONCOVALENT) -> pd.DataFrame:
    """Top pose by ONE chosen GNINA signal within each (pocket, compound, docking mode). Other signals are carried along unmodified.
    Covalent and noncovalent rows are never ranked together; pass mode=None to keep both lanes (each ranked within its own mode)."""
    d = res[res.processing_status == "docked"]
    if mode is not None: d = d[d.docking_mode == mode]
    # covalent lane: rank only poses where the modelled covalent bond is satisfied; invalid poses are never shown as a compound's best
    d = d[(d.docking_mode != MODE_COVALENT) | (d.covalent_geometry == "covalent_bond_modelled")]
    asc = sort_by == "gnina_minimized_affinity_kcal"
    d = d.sort_values(sort_by, ascending=asc).groupby(["pocket_id", "ligand_key", "docking_mode"], as_index=False).head(1)
    return d.sort_values(["docking_mode", "pocket_id", sort_by], ascending=[True, True, asc])
