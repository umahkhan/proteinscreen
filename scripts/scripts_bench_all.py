import sys, json
import pandas as pd
from pscreen.docking import DockConfig, MODE_COVALENT
from pscreen.validation import redock, covalent_bench as cb
P = lambda n,t,c,s: print(f"[{n}/{t}] {c}: {s}", flush=True)
ids = ["1HSG","1STP","3ERT","1M17","1EVE","1KE5","1HWI","1GPK","1JD0","1L2S","1HQ2","4DFR"]
print("=== REDOCK v2 (receptor = all chains within 5 A of ligand)", flush=True)
r = redock.run_redock(ids, "runs/redock_v2", DockConfig(exhaustiveness=8, num_modes=9, cpu=3, workers=2, timeout_s=1200), progress=P)
print(r.drop(columns=["reason"]).to_string()); print(json.dumps(redock.summarize(r), indent=1), flush=True)
E = cb.select_entries(24, seed=0, heavy=(10, 30))
cfg = DockConfig(exhaustiveness=8, num_modes=9, cpu=3, workers=2, timeout_s=900)
print("=== COVALENT form A (deposited graph), both modes", len(E), flush=True)
a = cb.run(E, "runs/covbench_A", cfg, progress=P, form="ccd")
print(cb.summarize(a).to_string(), flush=True)
print("=== COVALENT form B (pscreen adduct) on entries where it applies, covalent only", flush=True)
try:
    b = cb.run(E, "runs/covbench_B", cfg, progress=P, form="pscreen", modes=(MODE_COVALENT,), only_adducts=True)
    print(cb.summarize(b).to_string(), flush=True)
except Exception as ex:
    print("form B skipped:", ex, flush=True)
print("ALL DONE", flush=True)
