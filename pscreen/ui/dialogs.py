"""Compound details popup: everything that is too much for the table."""
from __future__ import annotations
import io
import pandas as pd
import streamlit as st
from rdkit import Chem

from ..docking import MODE_COVALENT
from ..viewer import pose_block
from .depict import compound_png, warhead_png
from .interactions import analyse
from .style import SCORE_LEGEND, badge
from .tables import admet_rows, hits_of
from .components import GEO_HELP, evidence_badges


@st.dialog("Compound details", width="large")
def compound_dialog(R: dict, item: dict, lane: str, in_top: bool, target: str):
    row = item["row"]; prow = R["prepared"][R["prepared"].ligand_key == row.ligand_key].iloc[0]; hits = hits_of(prow)
    m = tuple(sorted({a for h in hits for a in h["match_atoms"]})); rx = tuple(sorted({a for h in hits for a in h["reactive_atoms"]}))
    st.markdown(f"### {row.compound_id}  ·  site `{row.pocket_id}`")
    st.markdown("".join(badge(t, k) for t, k in evidence_badges(row, lane, in_top)), unsafe_allow_html=True)
    c1, c2 = st.columns([1.2, 1])
    with c1:
        st.image(compound_png(prow.ligand_sdf, m, rx, 420, 280), caption="Warhead atoms shaded; reactive atom in red")
        st.code(prow.canonical_smiles, language=None)
    with c2:
        if hits:
            st.markdown("**Warhead(s)**  — chemical annotation, not proof of reactivity")
            for h in hits:
                a, b = st.columns([1, 1.3])
                a.image(warhead_png(prow.ligand_sdf, tuple(h["match_atoms"]), tuple(h["reactive_atoms"]), 200, 140))
                b.markdown(f"**{h['warhead'].replace('_', ' ')}**  \n{h['mechanism']}  \nReversible: {'yes' if h['reversible'] else 'no'}  \nTypical residues: {', '.join(h['typical_residues'])}")
        else: st.info("No recognised warhead motif. Eligible for noncovalent ranking only.")
        if item["geometry"]:
            g = item["geometry"]; st.markdown(f"**Reactive atom → {target}:** {g['distance']} Å · *{g['label']}*  \n{GEO_HELP.get(g['label'], '')}")
    st.markdown("**Docking signals**  (separate evidence; never summed)")
    st.dataframe(pd.DataFrame([{"Signal": a, "Value": row[c], "Meaning": b} for a, b, c in SCORE_LEGEND]), hide_index=True, use_container_width=True,
                 column_config={"Value": st.column_config.NumberColumn(format="%.3f", width="small")})
    res = R["results"]; d = res[(res.processing_status == "docked") & (res.compound_id == row.compound_id) & (res.docking_mode == lane)]
    if lane == MODE_COVALENT: d = d[d.covalent_geometry == "covalent_bond_modelled"]
    d = d.sort_values("gnina_cnn_vs", ascending=False).drop_duplicates("pocket_id")
    if len(d) > 1:
        st.markdown("**This compound at each site**  (best pose per site)")
        st.dataframe(d[["pocket_id", "gnina_cnn_pose_score", "gnina_minimized_affinity_kcal", "gnina_cnn_affinity_pK", "reactive_to_target_distance"]].rename(columns={
            "pocket_id": "Site", "gnina_cnn_pose_score": "CNN pose", "gnina_minimized_affinity_kcal": "Vina-like", "gnina_cnn_affinity_pK": "CNN pK", "reactive_to_target_distance": "Reactive→target Å"}),
            hide_index=True, use_container_width=True)
    if lane == MODE_COVALENT: st.caption("Covalent-lane scores are from GNINA covalent docking of the adduct form and are not comparable with noncovalent scores.")
    ad = admet_rows(row, R.get("admet_meta"))
    if len(ad):
        st.markdown("**Predicted ADMET properties** (ADMET-AI · Chemprop-based ensemble)")
        st.dataframe(ad, hide_index=True, use_container_width=True, column_config={"Prediction": st.column_config.NumberColumn(format="%.3f")})
        st.caption("Classification outputs are probabilities, not calibrated risks. These models do not predict target binding, covalent reactivity, potency or selectivity.")
    pose = pose_block(item["pose_sdf"], item["pose_idx"]); df, s = analyse(item["pocket_receptor"], pose)
    st.markdown("**Pose quality and contacts** (geometric heuristics)")
    k = st.columns(5); k[0].metric("Residues in contact", s["residues_in_contact"]); k[1].metric("H-bond candidates", s["hbond_candidate_residues"])
    k[2].metric("Hydrophobic", s["hydrophobic_residues"]); k[3].metric("Clashes <2.5 Å", s["clashes_lt_2.5A"]); k[4].metric("Closest (Å)", s["closest_protein_ligand_A"])
    with st.expander(f"Residue contacts ({len(df)})"): st.dataframe(df, hide_index=True, use_container_width=True)
    with st.expander(f"All retained poses ({len(item['all_poses'])})"):
        st.dataframe(item["all_poses"][["pose_idx", "gnina_minimized_affinity_kcal", "gnina_cnn_pose_score", "gnina_cnn_affinity_pK", "gnina_cnn_vs", "reactive_to_target_distance"]], hide_index=True, use_container_width=True)
    buf = io.StringIO(); w = Chem.SDWriter(buf); mm = Chem.Mol(pose); mm.SetProp("_Name", f"{row.compound_id}|{row.pocket_id}"); w.write(mm); w.close()
    st.download_button("Download this pose (SDF)", buf.getvalue(), f"{row.compound_id}_{row.pocket_id}.sdf")
