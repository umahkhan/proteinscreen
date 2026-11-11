import json, os, time
from pathlib import Path
import importlib, pytest


@pytest.fixture
def hosted(tmp_path, monkeypatch):
    monkeypatch.setenv("PSCREEN_HOSTED", "1"); monkeypatch.setenv("PSCREEN_DATA", str(tmp_path)); monkeypatch.setenv("PSCREEN_RETENTION_H", "1")
    import pscreen.hosting as h
    yield importlib.reload(h)
    for k in ("PSCREEN_HOSTED", "PSCREEN_DATA", "PSCREEN_RETENTION_H"): monkeypatch.delenv(k, raising=False)
    importlib.reload(h)                                  # never leak hosted mode into other tests


def test_workspaces_are_isolated_and_seeded(hosted):
    a, b = hosted.new_token(), hosted.new_token()
    ra, rb = hosted.ensure_workspace(a), hosted.ensure_workspace(b)
    assert ra != rb and ra.parent.name == a
    (ra / "secret_run").mkdir(); (ra / "secret_run" / "run_config.json").write_text("{}")
    assert not (rb / "secret_run").exists()
    assert any(p.name.startswith("example_") for p in ra.iterdir())          # demo run present
    cfg = (ra / "example_WRN_Cys727" / "run_config.json").read_text()
    assert "@RUN@" not in cfg and str(ra) in cfg                             # placeholder resolved to the real folder


def test_bad_tokens_rejected(hosted):
    for t in (None, "", "../../etc", "x" * 32, "G" * 32): 
        assert not hosted.valid_token(t)
    with pytest.raises(ValueError): hosted.runs_root("../../etc")


def test_caps(hosted):
    ok = dict(target_residue="A:CYS1", target_residues=[])
    assert hosted.check_config(ok, 10, 1000, 1000, 0, 0) == []
    assert hosted.check_config(ok, 999, 1000, 1000, 0, 0)                     # too many compounds
    assert hosted.check_config(ok, 10, 1000, 1000, 1, 0)                      # already one active run
    assert hosted.check_config(ok, 10, 99 * 2**20, 1000, 0, 0)                # huge protein
    k = hosted.clamp_config(dict(exhaustiveness=64, cpu_per_job=16, workers=8, timeout_s=99999, max_pockets_per_tool=9))
    assert k["exhaustiveness"] <= 8 and k["cpu_per_job"] <= 2 and k["workers"] == 1 and k["timeout_s"] <= 600 and k["max_pockets_per_tool"] <= 2


def test_purge_old_removes_idle_only(hosted):
    old, new = hosted.new_token(), hosted.new_token()
    hosted.ensure_workspace(old); hosted.ensure_workspace(new)
    (hosted.workspace_dir(old) / ".active").write_text(str(time.time() - 5 * 3600))
    assert hosted.purge_old() == [old]
    assert hosted.workspace_dir(new).exists() and not hosted.workspace_dir(old).exists()


def test_slot_is_exclusive(hosted, monkeypatch):
    monkeypatch.setattr(hosted, "SLOTS", 1)
    s1 = hosted.Slot(); s1.acquire()
    s2 = hosted.Slot(); calls = []
    import threading
    t = threading.Thread(target=lambda: s2.acquire(on_wait=lambda: calls.append(1), poll_s=0.05)); t.start()
    time.sleep(0.3); assert calls and t.is_alive()                            # second run waits
    s1.release(); t.join(3); assert not t.is_alive()


def test_local_mode_unchanged(monkeypatch):
    monkeypatch.delenv("PSCREEN_HOSTED", raising=False)
    import pscreen.hosting as h; h = importlib.reload(h)
    assert h.runs_root(None) == Path("runs") and h.check_config({}, 10**6, 10**9, 10**9, 9, 99) == []


def test_many_pockets_of_one_model_keep_valid_colours():
    import re
    from pscreen.ui import viewer3d
    pk = [{"id": f"lacuna_{i}", "residues": ["A:ALA1"], "center": [0, 0, 0], "size": [20, 20, 20], "source": "predicted_cryptic"} for i in range(1, 12)]
    assert all(re.fullmatch(r"#[0-9a-f]{6}", e["color"]) for e in viewer3d._entries(pk, None))


def test_safe_name_blocks_path_tricks(hosted):
    for bad in ("../../etc/passwd", "..\\..\\x.pdb", "/abs/path/x.pdb", ".hidden", "a b$c;rm -rf.pdb", "", "..", "x" * 300 + ".pdb"):
        n = hosted.safe_name(bad)
        assert "/" not in n and "\\" not in n and not n.startswith(".") and len(n) <= 80 and n
    assert hosted.safe_name("my protein.pdb") == "my_protein.pdb" and hosted.safe_name("3MXF.pdb") == "3MXF.pdb"


def test_global_queue_cap(hosted, monkeypatch):
    monkeypatch.setattr(hosted, "MAX_QUEUE", 2)
    ok = dict(target_residue="", target_residues=[])
    for i in range(2):
        t = hosted.new_token(); r = hosted.ensure_workspace(t) / f"r{i}"; r.mkdir(); (r / "pid").write_text(str(os.getpid()))
    assert any("busy" in m for m in hosted.check_config(ok, 5, 1000, 1000, 0, 0))
