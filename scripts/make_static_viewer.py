"""Build the free static viewer (deploy/static/): a read-only 3D explorer + results table for the shipped WRN Cys727 example.
No server needed: Mol* is loaded from a CDN, the structure and pose are embedded in explorer.html.   python scripts/make_static_viewer.py"""
import html, shutil, sys, tempfile
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pscreen import hosting
from pscreen.runio import load_run
from pscreen.docking import MODE_COVALENT
from pscreen.ui import viewer3d
from pscreen.ui.components import build_scene

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "deploy" / "static"; OUT.mkdir(parents=True, exist_ok=True)
tmp = Path(tempfile.mkdtemp()) / "ex"; shutil.copytree(ROOT / "examples" / "example_WRN_Cys727", tmp); hosting.relocate(tmp)
R = load_run(tmp); target = R["docked_target"]
pk = {p.pocket_id: p for p in R["pockets"]}
site = pk["lacuna_3"]
pdict = lambda p: {"id": p.pocket_id, "residues": p.residues, "center": list(p.center), "size": list(p.size), "source": p.source,
                   "residue_scores": p.metrics.get("residue_scores"), "residue_score_label": p.metrics.get("residue_score_label", "")}
items = build_scene(R, ["VVD-133214"], site, MODE_COVALENT, target, {})
assert items, "no covalent pose found"
rec = site.receptor_path or R["receptor"]
panel = viewer3d.panel("WRN helicase · apo structure · Cys727", f"covalent pose of VVD-133214 at {site.pocket_id} (Lacuna conformer) · predicted, not evidence of binding", rec,
                       [pdict(site), pdict(pk["pocketminer_1"])], items, target, None, site.pocket_id)
(OUT / "explorer.html").write_text(viewer3d.explorer_html([panel]))

res = R["results"]; d = res[res.processing_status == "docked"]
cov = d[(d.docking_mode == MODE_COVALENT) & (d.covalent_geometry == "covalent_bond_modelled")]
top = (cov.sort_values("gnina_cnn_vs", ascending=False).drop_duplicates(["compound_id", "pocket_id"])
       [["compound_id", "pocket_id", "gnina_cnn_vs", "gnina_cnn_pose_score", "gnina_minimized_affinity_kcal", "reactive_to_target_distance"]])
top.columns = ["Compound", "Site", "CNN VS", "CNN pose", "Vina-like (kcal/mol)", "Warhead to Cys727 (Å)"]
rows = "".join("<tr>" + "".join(f"<td>{html.escape(f'{v:.2f}' if isinstance(v, float) else str(v))}</td>" for v in r) + "</tr>" for r in top.itertuples(index=False))
head = "".join(f"<th>{c}</th>" for c in top.columns)
(OUT / "index.html").write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>pscreen · pocket finding and compound prioritization</title>
<style>body{{margin:0;font:16px/1.5 system-ui,sans-serif;color:#1b1b1a;background:#fff}}main{{max-width:1100px;margin:auto;padding:28px 20px}}
h1{{margin:0 0 4px}}.sub{{color:#55544f;margin:0 0 18px}}.btn{{display:inline-block;background:#2a78d6;color:#fff;padding:10px 18px;border-radius:8px;text-decoration:none;font-weight:700}}
iframe{{width:100%;height:640px;border:1px solid #ddd;border-radius:10px;background:#14143a}}table{{border-collapse:collapse;width:100%;font-size:14px}}
th,td{{padding:6px 10px;border-bottom:1px solid #e5e5e0;text-align:left}}th{{background:#f6f6f4}}.note{{background:#fff7e8;border:1px solid #f0d9a8;border-radius:8px;padding:12px 14px;margin:18px 0}}
li{{margin:4px 0}}</style></head><body><main>
<h1>pscreen</h1><p class="sub">Decide which compounds to test first against a protein: find pockets (including cryptic ones), dock compounds (including covalent poses), rank, and explore in 3D.</p>
<ul><li>Upload a protein and a CSV of compounds (<code>id,smiles</code>).</li><li>Pockets from P2Rank, PocketMiner and Lacuna; GNINA docking with a covalent lane; ADMET-AI properties.</li>
<li>Everything runs in a visual workflow, repeatable, with a full export.</li></ul>
<p><b>Status:</b> this page is a read-only example. Running pscreen on your own protein and compounds is currently available only as a local install; a public release is planned.</p>
<h2>Example: WRN helicase, Cys727</h2>
<p>An apo (closed) structure, a Vividion-type covalent compound, and the pose pscreen modelled bonded to Cys727. Drag to rotate, scroll to zoom.</p>
<iframe src="explorer.html" title="3D explorer" loading="lazy"></iframe>
<h3>Covalent poses (best per compound and site)</h3><table><tr>{head}</tr>{rows}</table>
<div class="note"><b>Predictions only.</b> Pockets and docked poses are computational hypotheses. A modelled covalent bond is not evidence of engagement, potency or selectivity.
Scores are weak when the pocket is closed in the input structure (note the positive Vina-like value above). The example run used a 17 Å limit around Cys727 for covalent docking, which biases results toward that site.</div>
<p style="color:#77766f;font-size:13px">Tools: P2Rank, PocketMiner, Lacuna, GNINA, RDKit, ADMET-AI, Mol*. See the project page for citations and licences.</p>
</main></body></html>""")
print("wrote", OUT, [p.name for p in OUT.iterdir()]); print(top.to_string(index=False))
