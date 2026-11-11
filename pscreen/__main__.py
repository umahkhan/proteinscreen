"""Command line:  python -m pscreen <check|run|export|list> ...   (the same pipeline the web app drives)"""
from __future__ import annotations
import argparse, json, sys, time
from . import api


def _common(p):
    p.add_argument("--protein", required=True, help="PDB or mmCIF file"); p.add_argument("--compounds", required=True, help=".csv/.tsv/.smi/.sdf with SMILES")
    p.add_argument("--residue", action="append", default=[], metavar="A:CYS12", help="covalent target residue (repeatable)")
    p.add_argument("--known-ligand", default="", help="PDB/SDF/MOL2 ligand defining the site"); p.add_argument("--box", default="", metavar="X,Y,Z[,SX,SY,SZ]", help="manual search box")


def _box(s):
    v = [float(x) for x in s.split(",")] if s else []
    return (v[:3], v[3:6] or None) if v else ([], None)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pscreen"); sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="validate inputs only"); _common(c)
    r = sub.add_parser("run", help="validate, run the pipeline, write runs/<name>/outputs/"); _common(r)
    r.add_argument("--name", required=True); r.add_argument("--pockets", default="", help="comma list of p2rank,pocketminer,lacuna (default p2rank when no known ligand/box)")
    r.add_argument("--admet", action="store_true"); r.add_argument("--no-covalent", action="store_true"); r.add_argument("--fast", action="store_true", help="exhaustiveness 4, 3 poses")
    r.add_argument("--timeout", type=int, default=0, help="seconds per docking job"); r.add_argument("--root", default="runs")
    r.add_argument("--pocketminer-threshold", type=float, default=0.0); r.add_argument("--covalent-max-center", type=float, default=0.0)
    e = sub.add_parser("export", help="(re)write the output package of a run"); e.add_argument("rundir")
    sub.add_parser("list", help="list runs")
    argv = list(sys.argv[1:] if argv is None else argv)
    for i, t in enumerate(argv[:-1]):                       # "--box -1.7,14.4,16.9": a leading minus sign must not be read as an option
        if t == "--box": argv[i:i + 2] = [f"--box={argv[i + 1]}"]; break
    a = ap.parse_args(argv)
    if a.cmd == "list":
        from .runio import list_runs
        from .ui.workflow import run_state
        for d in list_runs(): print(f"{run_state(d)[0]:12s} {d.name}")
        return 0
    if a.cmd == "export": print(api.export_outputs(a.rundir)); return 0
    ctr, size = _box(a.box)
    pf = api.preflight(a.protein, a.compounds, a.residue, a.known_ligand, ctr)
    for w in pf["warnings"]: print("warning:", w)
    for e_ in pf["errors"]: print("ERROR:", e_)
    if a.cmd == "check" or not pf["ok"]: print(json.dumps(pf["info"])); return 0 if pf["ok"] else 2
    extra = {k: v for k, v in (("timeout_s", a.timeout), ("pocketminer_threshold", a.pocketminer_threshold), ("covalent_max_center_A", a.covalent_max_center)) if v}
    rd = api.create_run(a.name, a.protein, a.compounds, a.root, a.residue, a.known_ligand, ctr, size, [t for t in a.pockets.split(",") if t], a.admet, not a.no_covalent, a.fast, **extra)
    from .pipeline import run
    t0 = time.time(); print(f"running {rd} ...", flush=True)
    try: run(rd)
    except Exception as ex: print(f"FAILED after {time.time()-t0:.0f}s: {ex}"); return 1
    print(f"finished in {time.time()-t0:.0f}s -> {rd/'outputs'}"); return 0


if __name__ == "__main__": sys.exit(main())
