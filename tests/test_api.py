from pathlib import Path
from pscreen import api

REC = "tests/fixtures/6AFC_rec.pdb"


def test_preflight_rejects_bad_inputs(tmp_path):
    assert not api.preflight("nope.pdb", "nope.csv")["ok"]
    lib = tmp_path / "l.csv"; lib.write_text("id,name\na,b\n")
    assert "No SMILES column" in " ".join(api.preflight(REC, str(lib))["errors"])
    lib.write_text("id,smiles\na,C=CC(=O)N\nb,xyz(\n")
    r = api.preflight(REC, str(lib)); assert r["ok"] and any("cannot be parsed" in w for w in r["warnings"])
    assert not api.preflight(REC, str(lib), residues=["A:CYS9999"])["ok"]


def test_create_run_and_slug(tmp_path):
    lib = tmp_path / "l.csv"; lib.write_text("id,smiles\na,CCO\n")
    rd = api.create_run("My run / 1", REC, str(lib), root=tmp_path, fast=True, pockets=["p2rank"])
    assert rd.name == "My_run_1" and (rd / "run_config.json").exists()
    try: api.create_run("My run / 1", REC, str(lib), root=tmp_path); assert False
    except FileExistsError: pass
