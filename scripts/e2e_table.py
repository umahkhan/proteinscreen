"""Tabulate every runs/e2e_* run into docs/E2E_RESULTS.md (status, time, sites, top compound per lane, problems)."""
import json
from pathlib import Path
import pandas as pd
rows, detail = [], []
for d in sorted(Path("runs").glob("e2e_*")):
    pr = json.loads((d / "progress.json").read_text()) if (d / "progress.json").exists() else {}
    sm = json.loads((d / "outputs" / "summary.json").read_text()) if (d / "outputs" / "summary.json").exists() else {}
    rk = pd.read_csv(d / "outputs" / "ranked_compounds.csv") if (d / "outputs" / "ranked_compounds.csv").exists() else pd.DataFrame()
    st = pd.read_csv(d / "outputs" / "compound_status.csv") if (d / "outputs" / "compound_status.csv").exists() else pd.DataFrame()
    secs = round(pr.get("ended", 0) - pr.get("started", 0)) if pr.get("ended") else None
    lib = sm.get("library") or {}
    pk = "; ".join(f"{k}: {v}" for k, v in (sm.get("pocket_tools") or {}).items())
    docked = int((st.final_status == "docked").sum()) if len(st) else 0
    top = "; ".join(f"{m}: {g.sort_values('rank').iloc[0].compound_id} (VS {g.sort_values('rank').iloc[0].gnina_cnn_vs:.2f})" for m, g in rk.groupby("docking_mode")) if len(rk) else "—"
    rows.append({"run": d.name[4:], "result": "SUCCESS" if pr.get("state") == "finished" else f"FAILED ({pr.get('state')})", "seconds": secs, "sites": pk, "compounds in": lib.get("total"),
                 "docked": docked, "not docked": (lib.get("total") or 0) - docked, "top compound per lane": top})
    for _, r in rk.iterrows(): detail.append({"run": d.name[4:], "lane": r.docking_mode, "rank": r["rank"], "compound": r.compound_id, "site": r.pocket_id, "CNN VS": round(r.gnina_cnn_vs, 2),
                                           "kcal/mol": round(r.gnina_minimized_affinity_kcal, 2), "warhead": r.get("warhead_classes") if pd.notna(r.get("warhead_classes")) else "",
                                           "cov. bond Å": round(r.reactive_to_target_distance, 2) if pd.notna(r.get("reactive_to_target_distance")) else ""})
    for _, r in st[st.final_status != "docked"].iterrows(): detail.append({"run": d.name[4:], "lane": "—", "rank": "", "compound": r.compound_id, "site": "", "CNN VS": "", "kcal/mol": "", "warhead": "", "cov. bond Å": f"{r.final_status}: {r.prep_reason}"})
T, D = pd.DataFrame(rows), pd.DataFrame(detail)
out = ["# End-to-end test results", "", "Each row is one run of `python -m pscreen run --fast ...` (GNINA exhaustiveness 4, 3 poses, CPU, 6-core VM). Scores are docking predictions, not affinities.", "",
       T.to_markdown(index=False), "", "## Per-compound results", "", D.to_markdown(index=False), ""]
Path("docs/E2E_RESULTS.md").write_text("\n".join(out)); print("\n".join(out))
