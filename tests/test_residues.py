from pscreen.residues import validate_target, parse_spec, residue_table, near_pocket

REC = "tests/fixtures/6AFC_rec.pdb"      # DJ-1, chain A, Cys106


def test_parse_variants():
    assert parse_spec("A:CYS145")[:3] == ("A", 145, "CYS")
    assert parse_spec("cys 145")[:3] == (None, 145, "CYS")
    assert parse_spec("A:145")[:3] == ("A", 145, None)
    assert parse_spec("nonsense")[3] is not None and parse_spec("")[3] is not None


def test_valid_and_mismatch_and_missing():
    ok = validate_target(REC, "A:CYS106"); assert ok["ok"] and ok["nucleophile_atom"] == "SG" and ok["spec"] == "A:CYS106"
    bad = validate_target(REC, "A:CYS107"); assert not bad["ok"] and "is ALA" in bad["errors"][0] and "CYS residues" in bad["errors"][0]  # 107 is Ala; hints list real cysteines
    gone = validate_target(REC, "A:CYS9999"); assert not gone["ok"] and "does not exist" in gone["errors"][0]
    chain = validate_target(REC, "B:CYS106"); assert not chain["ok"] and "Chain 'B'" in chain["errors"][0]


def test_unsupported_residue_and_pocket_distance_warning():
    gly = validate_target(REC, "A:GLY74"); assert not gly["ok"] and "no nucleophilic" in gly["errors"][0]
    far = validate_target(REC, "A:CYS106", pocket_center=(100.0, 100.0, 100.0)); assert far["ok"] and far["warnings"] and far["distance_to_pocket"] > 15


def test_near_pocket_lists_reactive_residues_sorted():
    t = residue_table(REC); c = t[t.resnum == 106].iloc[0].nuc_xyz
    n = near_pocket(t, c, 12.0); assert n.iloc[0].resnum == 106 and list(n.dist_to_pocket) == sorted(n.dist_to_pocket)
