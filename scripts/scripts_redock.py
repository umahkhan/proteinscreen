import sys, json
from pscreen.docking import DockConfig
from pscreen.validation.redock import run_redock, summarize
ids = ["1HSG","1STP","3ERT","1M17","1EVE","1KE5","1HWI","1GPK","1JD0","1L2S","1HQ2","4DFR"]
res = run_redock(ids, "runs/redock_v1", DockConfig(exhaustiveness=8, num_modes=9, cpu=3, workers=2, timeout_s=1500),
                 progress=lambda n,t,c,s: print(f"[{n}/{t}] {c}: {s}", flush=True))
print(res.to_string()); print(json.dumps(summarize(res), indent=1))
