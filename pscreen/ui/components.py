"""Reusable UI pieces: scene building, the three-engine structure viewer, compound card, charts."""
from __future__ import annotations
import json
import altair as alt
import pandas as pd
from rdkit import Chem

from ..docking import MODE_COVALENT
from ..geometry import load_structure, residue_atoms, reactive_geometry
from ..viewer import pose_block
from .style import ORANGE, LANE_COLOR

SIGNAL_LABEL = {"gnina_minimized_affinity_kcal": "Vina-like score (kcal/mol, lower = better)", "gnina_cnn_pose_score": "CNN pose score (0–1)",
                "gnina_cnn_affinity_pK": "CNN affinity (pK)", "gnina_cnn_vs": "CNN VS (pose × affinity)"}
GEO_HELP = {"plausible": "reactive atom ≤ 4 Å from the nucleophile in a noncovalent pose", "marginal": "4–6 Å", "unlikely": "farther than 6 Å", "unmapped": "pose could not be mapped to the ligand",
            "not_assessed": "no target residue / no warhead", "covalent_bond_modelled": "covalent-lane pose with the modelled bond satisfied (≤ 2.3 Å)", "bond_not_formed": "bond not satisfied"}


# ------------------------------------------------------------------------------------------------------------- scene building
def build_scene(R: dict, compounds: list[str], pocket, lane: str, target: str, pose_choice: dict[str, int]) -> list[dict]:
    """One entry per compound with REAL pose coordinates (GNINA output), warhead atoms and reactive-atom geometry."""
    res, prepared, cmap = R["results"], R["prepared"], R["covalent_map"]
    docked = res[res.processing_status == "docked"]
    prec = pocket.receptor_path or R["receptor"]
    try: tgt = residue_atoms(load_structure(prec), target) if target else None
    except ValueError: tgt = None
    out = []
    for c in compounds:
        pp = docked[(docked.compound_id == c) & (docked.pocket_id == pocket.pocket_id) & (docked.docking_mode == lane)]
        if lane == MODE_COVALENT: pp = pp[pp.covalent_geometry == "covalent_bond_modelled"]
        pp = pp.sort_values("gnina_cnn_vs", ascending=False)
        if pp.empty: continue
        idx = int(pose_choice.get(c, int(pp.iloc[0].pose_idx)))
        row = pp[pp.pose_idx == idx].iloc[0] if (pp.pose_idx == idx).any() else pp.iloc[0]
        prow = prepared[prepared.ligand_key == row.ligand_key].iloc[0]
        pose = pose_block(row.pose_sdf, int(row.pose_idx))
        if lane == MODE_COVALENT and cmap is not None:
            cm = cmap[(cmap.ligand_key == row.ligand_key) & (cmap.warhead == row.geometry_warhead)].iloc[0]
            lig_sdf = cm.adduct_sdf; hits = [{"warhead": cm.warhead, "match_atoms": [int(cm.attach_idx)], "reactive_atoms": [int(cm.attach_idx)]}]
        else:
            lig_sdf = prow.ligand_sdf; hits = json.loads(prow.warheads) if isinstance(prow.get("warheads"), str) else []
        geom = None
        if tgt and hits:
            lig = Chem.RemoveHs(Chem.SDMolSupplier(lig_sdf, removeHs=False)[0])
            g = min((reactive_geometry(lig, pose, h["reactive_atoms"], tgt) for h in hits), key=lambda g: g.distance if g.distance is not None else 1e9)
            if g.distance is not None: geom = {"distance": g.distance, "reactive_xyz": g.reactive_xyz, "target_xyz": g.target_xyz, "label": g.label}
        out.append({"label": c, "pose_sdf": row.pose_sdf, "pose_idx": int(row.pose_idx), "ligand_sdf": lig_sdf, "warheads": hits, "geometry": geom, "row": row,
                    "pocket_receptor": prec, "all_poses": pp})
    return out


# ------------------------------------------------------------------------------------------------------------- viewer panel


# ------------------------------------------------------------------------------------------------------------- compound card
def has_wh(v) -> bool:
    return isinstance(v, str) and v.strip() not in ("", "nan", "None")


def evidence_badges(row, lane: str, in_top: bool) -> list[tuple[str, str]]:
    b = []
    if in_top: b.append(("top 10% by selected signal", "blue"))
    if has_wh(row.get("warhead_classes")): b.append((f"warhead motif: {row['warhead_classes']}", "orange"))
    g = row.get("covalent_geometry")
    if g in ("plausible", "covalent_bond_modelled"): b.append(("plausible covalent geometry", "orange"))
    elif g in ("marginal",): b.append(("marginal covalent geometry", "grey"))
    if any(str(c).startswith("exp_") and pd.notna(row.get(c)) for c in row.index): b.append(("has experimental data", "green"))
    if lane == MODE_COVALENT: b.append(("covalent-lane score (not comparable)", "orange"))
    return b or [("no flags", "grey")]


# ------------------------------------------------------------------------------------------------------------- charts
def _axis(title): return alt.Axis(title=title, labelColor="#5f5e5a", titleColor="#5f5e5a", gridColor="#eeeeea", domainColor="#d8d7d2")


def top_bar(best: pd.DataFrame, signal: str, lane: str, n: int = 15):
    asc = signal == "gnina_minimized_affinity_kcal"
    d = best.sort_values(signal, ascending=asc).head(n)
    return alt.Chart(d).mark_bar(color=LANE_COLOR[lane], size=12, cornerRadiusEnd=3).encode(
        y=alt.Y("compound_id:N", sort=alt.EncodingSortField(signal, order="ascending" if asc else "descending"), axis=alt.Axis(title=None, labelColor="#1b1b1a", labelLimit=200)),
        x=alt.X(f"{signal}:Q", axis=_axis(SIGNAL_LABEL[signal])),
        tooltip=["compound_id", "pocket_id", alt.Tooltip("gnina_minimized_affinity_kcal:Q", format=".2f"), alt.Tooltip("gnina_cnn_pose_score:Q", format=".3f"),
                 alt.Tooltip("gnina_cnn_affinity_pK:Q", format=".2f"), "warhead_classes"]).properties(height=26 * len(d) + 90, title=f"Top {len(d)} compounds")


def histogram(best: pd.DataFrame, signal: str, lane: str):
    return alt.Chart(best).mark_bar(color=LANE_COLOR[lane], opacity=0.9, cornerRadiusEnd=2).encode(
        x=alt.X(f"{signal}:Q", bin=alt.Bin(maxbins=24), axis=_axis(SIGNAL_LABEL[signal])), y=alt.Y("count()", axis=alt.Axis(title="Compounds", tickMinStep=1, labelColor="#5f5e5a", titleColor="#5f5e5a", gridColor="#eeeeea", domainColor="#d8d7d2")),
        tooltip=[alt.Tooltip("count()", title="compounds")]).properties(height=290, title="Distribution of the selected signal")


def size_scatter(best: pd.DataFrame, signal: str, lane: str):
    pts = alt.Chart(best).mark_circle(size=70, color=LANE_COLOR[lane], opacity=0.75, stroke="#ffffff", strokeWidth=1.5).encode(
        x=alt.X("mol_weight:Q", axis=_axis("Molecular weight")), y=alt.Y(f"{signal}:Q", axis=_axis(SIGNAL_LABEL[signal]), scale=alt.Scale(zero=False)),
        tooltip=["compound_id", alt.Tooltip("mol_weight:Q", format=".1f"), alt.Tooltip(f"{signal}:Q", format=".3f"), "warhead_classes"])
    return (pts + pts.transform_regression("mol_weight", signal).mark_line(color="#8a8a85", strokeDash=[4, 3])).properties(height=330, title="Size-bias check: score vs molecular weight")


def warhead_bar(res: pd.DataFrame):
    s = res[res.processing_status == "docked"].drop_duplicates("ligand_key").warhead_classes.fillna("").str.split(";").explode()
    s = s[s != ""].value_counts().rename_axis("class").reset_index(name="compounds")
    if s.empty: return None
    return alt.Chart(s).mark_bar(color=ORANGE, size=12, cornerRadiusEnd=3).encode(
        y=alt.Y("class:N", sort="-x", axis=alt.Axis(title=None, labelLimit=220)), x=alt.X("compounds:Q", axis=alt.Axis(title="Compounds", tickMinStep=1, labelColor="#5f5e5a", titleColor="#5f5e5a", gridColor="#eeeeea", domainColor="#d8d7d2")),
        tooltip=["class", "compounds"]).properties(height=28 * len(s) + 100, title="Recognised warhead classes")
