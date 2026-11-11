import pandas as pd, pytest
from pathlib import Path
from pscreen.results import best_pose_per_compound
from pscreen.pipeline import prepare_receptor


def _res():
    base = dict(processing_status="docked", pocket_id="P1", gnina_minimized_affinity_kcal=-5.0, gnina_cnn_pose_score=0.5, gnina_cnn_affinity_pK=5.0)
    rows = [  # compound A: best covalent score has the bond NOT formed -> must be ignored in the covalent lane
        dict(base, ligand_key="A", docking_mode="covalent", gnina_cnn_vs=9.0, covalent_geometry="bond_not_formed"),
        dict(base, ligand_key="A", docking_mode="covalent", gnina_cnn_vs=2.0, covalent_geometry="covalent_bond_modelled"),
        # compound B: covalent poses all invalid -> must not appear in covalent lane at all
        dict(base, ligand_key="B", docking_mode="covalent", gnina_cnn_vs=8.0, covalent_geometry="bond_not_formed"),
        dict(base, ligand_key="A", docking_mode="noncovalent", gnina_cnn_vs=1.0, covalent_geometry="not_assessed"),
        dict(base, ligand_key="A", docking_mode="noncovalent", gnina_cnn_vs=3.0, covalent_geometry="not_assessed"),
    ]
    return pd.DataFrame(rows)


def test_covalent_lane_only_ranks_bond_satisfied_poses():
    c = best_pose_per_compound(_res(), "gnina_cnn_vs", "covalent")
    assert list(c.ligand_key) == ["A"] and c.iloc[0].gnina_cnn_vs == 2.0


def test_lanes_are_not_mixed_and_noncovalent_unaffected():
    n = best_pose_per_compound(_res(), "gnina_cnn_vs", "noncovalent")
    assert set(n.docking_mode) == {"noncovalent"} and n.iloc[0].gnina_cnn_vs == 3.0
    both = best_pose_per_compound(_res(), "gnina_cnn_vs", None)
    assert set(both.docking_mode) == {"covalent", "noncovalent"} and len(both) == 2   # one best pose per (compound, lane); B dropped from covalent


def test_vina_like_signal_sorts_ascending():
    d = pd.DataFrame([dict(processing_status="docked", pocket_id="P", ligand_key="X", docking_mode="noncovalent", pose_idx=i, gnina_minimized_affinity_kcal=v,
                           gnina_cnn_vs=0.0, covalent_geometry="x") for i, v in enumerate([-3.0, -9.0, -5.0])])
    assert best_pose_per_compound(d, "gnina_minimized_affinity_kcal", "noncovalent").iloc[0].pose_idx == 1


def test_prepare_receptor_rejects_non_structure(tmp_path):
    bad = tmp_path / "x.pdb"; bad.write_text("not a structure\n")
    with pytest.raises(ValueError):
        prepare_receptor(bad, tmp_path / "out.pdb")


def test_prepare_receptor_strips_hetero_and_waters(tmp_path):
    src = Path("data/smoke/rec.pdb")          # protein only; add a ligand + water to be stripped
    txt = src.read_text().replace("END\n", "") + "HETATM 9000  O   HOH A 300      10.000  10.000  10.000  1.00  0.00           O\nHETATM 9001  C1  MK1 A 400      11.000  11.000  11.000  1.00  0.00           C\nEND\n"
    p = tmp_path / "in.pdb"; p.write_text(txt)
    info = prepare_receptor(p, tmp_path / "out.pdb")
    assert info["n_residues"] == 198 and "HETATM" not in (tmp_path / "out.pdb").read_text()
