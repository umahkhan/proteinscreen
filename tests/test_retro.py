import numpy as np, pandas as pd
from pscreen.validation.retro import evaluate


def test_perfect_and_inverted_signals_and_orientation():
    n = 40
    ids = [f"c{i}" for i in range(n)]
    y = np.array([1] * 10 + [0] * 30)
    # vina-like: more negative = better, so actives get lower values -> should score AUC 1.0 after negation
    best = pd.DataFrame({"compound_id": ids, "gnina_minimized_affinity_kcal": np.where(y == 1, -10.0, -5.0) + np.arange(n) * 1e-3,
                         "gnina_cnn_pose_score": np.where(y == 1, 0.2, 0.9), "mol_weight": 300.0 + np.arange(n)})
    r = evaluate(best, pd.DataFrame({"compound_id": ids, "active": y}), "active").set_index("signal")
    assert r.loc["gnina_minimized_affinity_kcal", "roc_auc"] == 1.0          # orientation handled
    assert r.loc["gnina_cnn_pose_score", "roc_auc"] == 0.0                   # inverted signal is reported as such, not flipped
    assert "baseline_mol_weight" in r.index


def test_continuous_labels_use_spearman_and_single_class_reported():
    ids = [f"c{i}" for i in range(20)]
    best = pd.DataFrame({"compound_id": ids, "gnina_cnn_affinity_pK": np.arange(20.0), "mol_weight": 1.0 * np.arange(20)})
    r = evaluate(best, pd.DataFrame({"compound_id": ids, "pic50": np.arange(20.0)}), "pic50", "continuous").set_index("signal")
    assert abs(r.loc["gnina_cnn_affinity_pK", "spearman"] - 1.0) < 1e-9
    one = evaluate(best, pd.DataFrame({"compound_id": ids, "a": [1] * 20}), "a").set_index("signal")
    assert one.loc["gnina_cnn_affinity_pK", "note"] == "only one class present"
