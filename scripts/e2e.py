"""End-to-end regression: small, diverse proteins x small SMILES sets, each through the product CLI (python -m pscreen run).
Usage: python scripts/e2e.py [case ...]   (no args = all, stops at the first failure so it can be fixed before moving on)
Results: runs/e2e_<case>/outputs/ and docs/E2E_RESULTS.md (table of every case run so far)."""
import json, subprocess, sys, time
from pathlib import Path
D = "data/e2e/"
CASES = {   # case: (description, args)
    "hivpr":        ("HIV-1 protease (1HSG, aspartyl protease), site from crystal ligand MK1, non-covalent", ["--protein", D+"1HSG.pdb", "--compounds", D+"hivpr.csv", "--known-ligand", D+"hivpr_MK1.pdb"]),
    "brd4":         ("BRD4 BD1 (3MXF, bromodomain), P2Rank + PocketMiner, ADMET, duplicate SMILES", ["--protein", D+"3MXF.pdb", "--compounds", D+"brd4.csv", "--pockets", "p2rank,pocketminer", "--admet"]),
    "streptavidin": ("Streptavidin (1STP, beta-barrel), P2Rank, one invalid SMILES", ["--protein", D+"1STP.pdb", "--compounds", D+"streptavidin.csv", "--pockets", "p2rank"]),
    "trypsin":      ("Trypsin (3PTB, serine protease), manual search box at benzamidine site", ["--protein", D+"3PTB.pdb", "--compounds", D+"trypsin.csv", "--box", "-1.76,14.46,16.92,18,18,18"]),
    "caii":         ("Carbonic anhydrase II (3HS4, Zn enzyme), manual box at acetazolamide site", ["--protein", D+"3HS4.pdb", "--compounds", D+"caii.csv", "--box", "-4.22,-4.08,17.04,18,18,18"]),
    "t4l":          ("T4 lysozyme L99A (181L, buried hydrophobic cavity), Lacuna ensemble pockets", ["--protein", D+"181L.pdb", "--compounds", D+"t4l.csv", "--pockets", "lacuna"]),
    "kras":         ("KRAS G12C (4LDJ chain A, GDP-bound), P2Rank + covalent docking at Cys12", ["--protein", D+"4LDJ_A.pdb", "--compounds", D+"kras_g12c.csv", "--pockets", "p2rank", "--residue", "A:CYS12",
                                                                                                   "--covalent-max-center", "14", "--timeout", "3000"]),
}


def run_case(k):
    desc, args = CASES[k]; name = f"e2e_{k}"
    t0 = time.time()
    p = subprocess.run([sys.executable, "-m", "pscreen", "run", "--name", name, "--fast", *args], capture_output=True, text=True)
    return {"case": k, "desc": desc, "rc": p.returncode, "seconds": round(time.time() - t0), "tail": (p.stdout + p.stderr).strip().splitlines()[-4:]}


if __name__ == "__main__":
    for k in (sys.argv[1:] or CASES):
        r = run_case(k); print(json.dumps(r, indent=1), flush=True)
        if r["rc"]: sys.exit(1)
