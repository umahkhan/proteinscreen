import json, sys
from pscreen.validation import cryptic
# 1) choose PocketMiner threshold on VAL only (top-1 pocketminer recovery; ties -> lowest threshold closest to default 0.7)
best = None
for thr in (0.5, 0.6, 0.7, 0.8):
    r = cryptic.run("val", f"runs/cryptic/val_thr{thr}", threshold=thr, tools=("pocketminer",))
    s = cryptic.summarize(r, thr); print("VAL thr", thr, json.dumps(s), flush=True)
    score = (s.get("pocketminer_top3_recovery", 0), -abs(thr - 0.7))
    if best is None or score > best[0]: best = (score, thr)
thr = best[1]; print("selected threshold (VAL top-3 recovery):", thr, flush=True)
r = cryptic.run("test", f"runs/cryptic/test_thr{thr}", threshold=thr, tools=("p2rank", "pocketminer", "lacuna"), progress=lambda i,s: print(i, s, flush=True))
s = cryptic.summarize(r, thr); print("TEST", json.dumps(s, indent=1)); json.dump({"threshold": thr, "test": s}, open("runs/cryptic/summary.json", "w"), indent=1)
print(r["status"][r["status"].status != "ok"].to_string())
