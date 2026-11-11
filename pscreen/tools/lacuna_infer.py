"""Runs in the isolated Lacuna env (lacuna-pockets, MIT). Mirrors `lacuna discover` defaults
(NMA backend, alpha detector, learned ranker) and additionally records, per pocket cluster, the conformer in which
that pocket is largest ("most open") and writes that conformer as its own PDB for docking.
Usage: lacuna_infer.py <in.pdb> <outdir> <n_conformers> <top> <seed>"""
import json, sys
from pathlib import Path
from lacuna import load_structure, detect_pockets, cluster_pockets
from lacuna.ensemble.nma_backend import NMABackend
from lacuna.io.structure import coords_array
from lacuna.io.writers import write_ensemble_pdb

inp, out, n, top, seed = sys.argv[1], Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
out.mkdir(parents=True, exist_ok=True)
structure = load_structure(inp)
coord_sets = NMABackend(seed=seed).generate(inp, n)
allc = [coords_array(structure)] + list(coord_sets)
pl = []
for ci, c in enumerate(allc):
    ps = detect_pockets(c, structure)
    for p in ps:
        p.conformer_idx = ci
    pl.append(ps)
clusters = cluster_pockets(pl, n_conformers=len(allc), rank_by="learned")[:top]
write_ensemble_pdb(structure, allc, out / "ensemble.pdb")

# split models so single conformers can be docked into
models, cur = {}, []
for line in open(out / "ensemble.pdb"):
    if line.startswith("MODEL"):
        cur = []; idx = int(line.split()[1]) - 1
    elif line.startswith("ENDMDL"):
        models[idx] = cur
    elif line.startswith(("ATOM", "HETATM", "TER")):
        cur.append(line)

report = []
for c in clusters:
    d = c.to_dict()
    best = max(c.member_pockets, key=lambda p: p.volume_a3)
    ci = int(best.conformer_idx)
    pdb = out / f"conformer_{ci}.pdb"
    if not pdb.exists():
        pdb.write_text("".join(models[ci]) + "END\n")
    d.update(volume_by_conformer={int(m.conformer_idx): round(float(m.volume_a3), 1) for m in c.member_pockets},
             open_conformer_idx=ci, open_conformer_pocket_volume_a3=float(best.volume_a3),
             open_conformer_pocket_centroid=[float(x) for x in best.centroid], open_conformer_pdb=str(pdb),
             input_conformer_pdb=str(inp))
    report.append(d)
json.dump({"seed": seed, "n_conformers": len(allc), "backend": "nma", "detector": "alpha", "ranked_by": "learned", "pockets": report},
          open(out / "report.json", "w"), indent=1, default=float)
