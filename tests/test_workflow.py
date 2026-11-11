import json, time
from pathlib import Path
from pscreen.pipeline import Progress, RunConfig, plan
from pscreen.ui import workflow as wf


def N(*types, **params):
    return [{"type": t, "id": f"n{i}", "x": 0, "y": 0, "params": params.get(t, {})} for i, t in enumerate(types)]


PATHS = {"receptor": "r.pdb", "library": "l.csv"}


def test_p2rank_only_no_admet():
    kw, err, _ = wf.to_config(N("protein", "library", "p2rank", "gnina", "warhead", "results"), "x_t1", PATHS, None)
    assert not err and kw["predictors"] == ["p2rank"] and not kw["chemprop"] and not kw["covalent"]


def test_two_pocket_models_with_admet_flag_follows_node(monkeypatch):
    monkeypatch.setattr(wf.chem_props, "available", lambda: True)
    monkeypatch.setattr(wf.pocket_tools, "PM_PY", Path(__file__))
    kw, err, _ = wf.to_config(N("p2rank", "pocketminer", "gnina", "admet", pocketminer={"threshold": 0.6, "max_pockets": 2}), "x_t2", PATHS, None)
    assert not err and kw["predictors"] == ["p2rank", "pocketminer"] and kw["chemprop"] and kw["pocketminer_threshold"] == 0.6


def test_covalent_residues_multiple_and_required():
    kw, err, _ = wf.to_config(N("p2rank", "gnina", "target", "gnina_cov", target={"residues": "A:CYS727, A:CYS693"}), "x_t3", PATHS, None)
    assert kw["covalent"] and kw["target_residue"] == "A:CYS727" and kw["target_residues"] == ["A:CYS693"]
    _, err, _ = wf.to_config(N("p2rank", "gnina", "target", "gnina_cov", target={"residues": ""}), "x_t4", PATHS, None)
    assert any("target residue" in e for e in err)


def test_validation_messages():
    _, err, _ = wf.to_config(N("gnina"), "x_t5", {}, None)
    assert len(err) >= 3                                       # protein, library, binding site
    _, err, _ = wf.to_config(N("p2rank"), "x_t6", PATHS, None)
    assert any("GNINA" in e for e in err)


def test_runconfig_residues_and_plan():
    c = RunConfig(name="a", receptor="r", library="l", target_residue="A:CYS1", target_residues=["A:CYS2", "A:CYS1"], predictors=["p2rank", "lacuna"], predict_pockets=True, covalent=True, chemprop=True)
    assert c.residues() == ["A:CYS1", "A:CYS2"]
    assert plan(c) == ["inputs", "prepare", "pockets", "pocket:p2rank", "pocket:lacuna", "docking", "covalent_docking", "properties", "results"]


def test_progress_times_each_stage(tmp_path):
    p = Progress(tmp_path, ["a", "b"]); p.set("a", msg="x"); time.sleep(0.05); p.set("b", done=1, total=2, msg="y"); time.sleep(0.05); p.finish("finished")
    d = json.loads((tmp_path / "progress.json").read_text()); a, b = d["stages"]["a"], d["stages"]["b"]
    assert a["end"] - a["start"] >= 0.05 and b["end"] >= b["start"] and a["status"] == b["status"] == "done" and d["ended"] >= d["started"]
