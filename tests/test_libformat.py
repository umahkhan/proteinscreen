from pscreen.libformat import TEMPLATE_CSV, validate_library


def test_template_is_valid(tmp_path):
    f = tmp_path / "t.csv"; f.write_text(TEMPLATE_CSV)
    r = validate_library(f); assert r.ok and r.n_ok == 3 and not r.warnings


def test_row_problems_are_reported_not_dropped(tmp_path):
    f = tmp_path / "m.csv"; f.write_text("id,smiles\na,CCO\na,CCN\nb,notasmiles\n,CCC\n")
    r = validate_library(f)
    assert r.n_total == 4 and r.n_ok == 2
    s = dict(zip(r.table["smiles"], r.table["status"]))
    assert "duplicate" in s["CCN"] and "SMILES" in s["notasmiles"]
    assert r.table.loc[r.table.smiles == "CCC", "id"].iloc[0] == "record_4"      # blank id numbered, not dropped
    assert r.standard_csv().splitlines()[0] == "id,smiles"


def test_file_level_errors_are_plain_language(tmp_path):
    f = tmp_path / "n.csv"; f.write_text("name,structure\nx,CCO\n"); r = validate_library(f)
    assert not r.ok and "smiles" in r.errors[0].lower() and "template" in r.errors[0]
    g = tmp_path / "x.pdf"; g.write_text("x"); assert "not supported" in validate_library(g).errors[0]
    h = tmp_path / "e.csv"; h.write_text("id,smiles\n"); assert not validate_library(h).ok
