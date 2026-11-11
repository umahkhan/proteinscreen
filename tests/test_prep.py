from pscreen.prep import prepare_library


def test_mini_library(tmp_path):
    df = prepare_library("data/examples/mini.csv", tmp_path).set_index("compound_id")
    assert len(df) == 9                                  # nothing dropped
    assert df.loc["aspirin_dup", "status"] == "duplicate" and df.loc["aspirin_dup", "duplicate_of"] == "aspirin"
    assert df.loc["bad_valence", "status"] == "invalid" and df.loc["garbage", "status"] == "invalid"
    assert df.loc["empty", "status"] == "invalid"
    assert "multiple fragments" in df.loc["salt_case", "warnings"]
    assert df.loc["acrylamide_aniline", "warhead_classes"] == "acrylamide"
    assert df.loc["chloroacetamide_1", "warhead_classes"] == "chloroacetamide"
    assert df.loc["benzamide", "status"] == "ok" and not df.loc["benzamide", "warhead_classes"]  # still dockable
    assert (tmp_path / "ligands" / "aspirin.sdf").exists()
    assert not (tmp_path / "ligands" / "aspirin_dup.sdf").exists()
