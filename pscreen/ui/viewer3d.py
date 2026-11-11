"""Structure explorer: a self-contained, full-height Mol* (molstar.org) page for inspecting predicted pockets and docked compounds together.

Visual language follows P2Rank's PyMOL output: dark gradient background, blue-white solvent-accessible surface (SAS), each predicted pocket a
coloured surface patch. Compounds are drawn as sticks where they docked and are identified by colour (key panel; hover for the name).

Highlighting (independent toggles, combinable):
  Ligand          translucent shell + thicker sticks around the pose(s)
  Pocket          solid coloured surface patch + thicker residue sticks and residue labels
  Target residue  yellow shell, thick sticks, label
While anything is highlighted the rest of the surface becomes faint and see-through, so the highlighted part is always visible, even when buried.
Residues can be added dynamically: click one in 3D, or tick residues listed near the visible pockets (radius / type filters).
All coordinates come from real files / GNINA output; nothing is illustrative.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from ..geometry import map_ligand_to_pose, load_structure
from ..viewer import pose_block

POSE_COLORS = ["cyan", "magenta", "green", "orange", "purple", "white"]       # 3Dmol '<name>Carbon' schemes
MODEL_COLORS = {"lacuna": "#3d7ff0", "p2rank": "#f28c28", "pocketminer": "#1fb872"}   # one hue per pocket-prediction approach
POCKET_COLORS = ["#5c93e6", "#e6935c", "#5ce693", "#bd5ce6", "#e65c93", "#e6d85c"]   # P2Rank's pocket palette
CDN = "https://cdn.jsdelivr.net/npm/molstar@5.13.1/build/viewer/molstar.js"
CSS = "https://cdn.jsdelivr.net/npm/molstar@5.13.1/build/viewer/molstar.css"
RES_CLASS = {"hydrophobic": ("#e0a030", ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP"]), "polar": ("#2fb5a0", ["SER", "THR", "ASN", "GLN", "TYR"]),
             "positive": ("#3b6fd6", ["LYS", "ARG", "HIS"]), "negative": ("#d94a4a", ["ASP", "GLU"]), "Gly / Pro": ("#9b8fd9", ["GLY", "PRO"]), "Cys": ("#e87ba4", ["CYS"])}
_CLS = {r: (c, n) for n, (c, rs) in RES_CLASS.items() for r in rs}
RAMP = ["#dbe8f8", "#a8c7ee", "#6fa2e0", "#2a78d6", "#1a4f93"]


def _resname(spec: str) -> str:
    import re
    return re.sub(r"[^A-Z]", "", spec.split(":")[1].upper())[:3]


def _ramp(v: float, lo: float, hi: float) -> str:
    t = 0.0 if hi <= lo else min(1.0, max(0.0, (v - lo) / (hi - lo))); return RAMP[min(len(RAMP) - 1, int(t * len(RAMP)))]


def polar_pairs(receptor_path: str, pose: Chem.Mol, cutoff: float = 3.5) -> list[list]:
    """N/O…N/O heavy-atom pairs within `cutoff` A (H-bond candidates; geometric heuristic). [[x1,y1,z1,x2,y2,z2,d,label], ...]"""
    P, xyz = [], []
    for model in load_structure(receptor_path):
        for ch in model:
            for r in ch:
                if r.id[0].strip() == "":
                    for a in r:
                        if a.element in ("N", "O"): P.append(f"{ch.id}:{r.get_resname()}{r.id[1]}:{a.get_id()}"); xyz.append(a.coord.astype(float))
        break
    if not P: return []
    tree = cKDTree(np.array(xyz)); conf = pose.GetConformer().GetPositions(); out = []
    for a in pose.GetAtoms():
        if a.GetSymbol() in ("N", "O"):
            for j in tree.query_ball_point(conf[a.GetIdx()], cutoff):
                d = float(np.linalg.norm(conf[a.GetIdx()] - xyz[j]))
                out.append([*map(float, conf[a.GetIdx()]), *map(float, xyz[j]), round(d, 2), P[j]])
    return out


def pose_entry(item: dict, receptor_path: str, color: str) -> dict:
    """item = one scene entry (components.build_scene) plus optional 'site'. Adds warhead/reactive atom coordinates, distance line, polar contacts."""
    pose = pose_block(item["pose_sdf"], int(item["pose_idx"]))
    heavy = Chem.RemoveHs(Chem.Mol(pose), sanitize=False)
    wh, rx = [], []
    if item.get("ligand_sdf") and item.get("warheads"):
        lig = Chem.RemoveHs(Chem.SDMolSupplier(item["ligand_sdf"], removeHs=False)[0])
        mp = map_ligand_to_pose(lig, pose); xyz = pose.GetConformer().GetPositions()
        if mp:
            for h in item["warheads"]:
                wh += [[float(v) for v in xyz[mp[a]]] for a in h["match_atoms"] if a in mp]
                rx += [[float(v) for v in xyz[mp[a]]] for a in h["reactive_atoms"] if a in mp]
    g = item.get("geometry")
    same_frame = item.get("pocket_receptor") in (None, receptor_path)      # polar-contact lines only where the pose was docked on the displayed receptor
    return {"label": item["label"], "site": item.get("site"), "data": Chem.MolToMolBlock(heavy, kekulize=False), "format": "sdf", "color": color, "warhead_xyz": wh, "reactive_xyz": rx,
            "distance": ({"a": list(g["reactive_xyz"]), "b": list(g["target_xyz"]), "text": f"{g['distance']:.1f} Å"} if g and g.get("reactive_xyz") else None),
            "polar": polar_pairs(receptor_path, pose) if same_frame else []}


def _model_key(pid: str) -> str:
    return next((k for k in MODEL_COLORS if pid.startswith(k)), "")


TINTS = [0.0, 0.28, 0.56, 0.68, 0.76, 0.82, 0.86, 0.90]       # lighter tint per extra pocket of one model (1st..3rd unchanged); capped so a colour never exceeds white


def _mix(hexc: str, t: float) -> str:
    t = min(1.0, max(0.0, t)); c = [int(hexc[j:j + 2], 16) for j in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (255 - x) * t):02x}" for x in c)


def _pocket_entry(p: dict, i: int, active: str | None, color: str | None = None) -> dict:
    scores = p.get("residue_scores") or {}
    vals = [scores[r] for r in p["residues"] if r in scores]
    lo, hi = (min(vals), max(vals)) if vals else (0, 1)
    return {**p, "model": _model_key(p["id"]), "color": color or POCKET_COLORS[i % len(POCKET_COLORS)], "active": (p["id"] == active),
            "res_class": {r: _CLS.get(_resname(r), ("#c9a000", "other"))[0] for r in p["residues"]},
            "res_score": {r: _ramp(scores[r], lo, hi) for r in p["residues"] if r in scores}, "scores": {r: scores[r] for r in p["residues"] if r in scores},
            "score_label": p.get("residue_score_label", ""), "score_range": [lo, hi]}


def _entries(pockets: list[dict], active: str | None) -> list[dict]:
    """Colour by prediction approach (Lacuna blue, P2Rank orange, PocketMiner green); extra pockets of one approach are lighter tints."""
    seen: dict[str, int] = {}; out = []
    for i, p in enumerate(pockets):
        k = _model_key(p["id"]); n = seen.get(k, 0); seen[k] = n + 1
        out.append(_pocket_entry(p, i, active, _mix(MODEL_COLORS[k], TINTS[min(n, len(TINTS) - 1)]) if k else None))
    return out


def panel(title: str, subtitle: str, receptor_path: str, pockets: list[dict], items: list[dict], target: str | None, reference_path: str | None,
          active_pocket: str | None = None, picked: list[str] | None = None, ensemble_path: str | None = None, ens_info: dict | None = None) -> dict:
    """pockets: ALL predicted/known sites to map [{id, residues, center, size, source?, residue_scores?}]; items: scene entries (poses) to overlay."""
    ref = None
    if reference_path and Path(reference_path).exists():
        ref = {"data": Path(reference_path).read_text(), "format": "pdb"}
    return {"title": title, "subtitle": subtitle, "receptor": Path(receptor_path).read_text(), "target": target or None, "reference": ref, "picked": picked or [],
            "pockets": _entries(pockets, active_pocket), "ensemble": Path(ensemble_path).read_text() if ensemble_path and Path(ensemble_path).exists() else None, "ens_info": ens_info,
            "poses": [pose_entry(it, receptor_path, POSE_COLORS[it.get("color_idx", i) % len(POSE_COLORS)]) for i, it in enumerate(items)]}


def explorer_html(panels: list[dict], dark: bool = True) -> str:
    """Full-height page (use components.html(..., height=H)); toolbar, viewport(s) and key panel always fit inside the iframe."""
    return TEMPLATE.replace("__SPEC__", json.dumps({"panels": panels, "dark": dark})).replace("__CDN__", CDN).replace("__CSS__", CSS)


TEMPLATE = r"""
<meta charset="utf-8"><link rel="stylesheet" href="__CSS__"/>
<style>
 :root{--ink:#1b1b1a;--mut:#6b6a66;--line:#e2e1dc;--soft:#f6f6f3;--acc:#2a78d6}
 *{box-sizing:border-box} html,body{height:100%;margin:0;font:13px/1.35 -apple-system,"Segoe UI",Inter,Helvetica,Arial,sans-serif;color:var(--ink);background:transparent}
 #root{position:relative;display:flex;flex-direction:column;height:100vh;border:1px solid var(--line);border-radius:12px;overflow:hidden;background:#fff;box-shadow:0 1px 3px rgba(0,0,0,.06)}
 #bar{display:flex;align-items:center;gap:10px;padding:8px 10px;border-bottom:1px solid var(--line);background:#fff;flex-wrap:wrap;position:relative;z-index:30}
 .grp{display:flex;align-items:center;gap:8px}.grp+.grp{padding-left:12px;border-left:1px solid var(--line)}
 .gl{font-size:11px;font-weight:800;letter-spacing:.07em;text-transform:uppercase;color:var(--ink)}
 button{font:inherit;color:inherit;cursor:pointer}
 .tool{padding:5px 11px;border:1.5px solid #8b93a1;border-radius:8px;background:#fff;font-weight:700;color:#111}.tool:hover{background:var(--soft)}.tool.open{background:#eef4fc;border-color:#9cbfe8}
 .seg{display:flex;border:1px solid var(--line);border-radius:9px;overflow:hidden}
 .seg button{border:0;border-right:1.5px solid #8b93a1;background:#fff;padding:6px 13px;display:flex;align-items:center;gap:7px;font-weight:700;color:#222}
 .seg button:last-child{border-right:0}.seg button:hover{background:var(--soft)}
 .seg i{width:10px;height:10px;border-radius:50%;background:var(--c);display:inline-block;opacity:.55}
 .seg button.on{background:var(--c);color:#fff;font-weight:700}.seg button.on i{background:#fff;opacity:1}
 .pop{display:none;position:absolute;top:100%;margin-top:6px;z-index:50;min-width:260px;max-width:340px;padding:12px;background:#fff;border:1px solid var(--line);border-radius:10px;box-shadow:0 8px 24px rgba(0,0,0,.14)}
 .pop.open{display:block}.pop h4{margin:0 0 8px;font-size:10px;letter-spacing:.07em;text-transform:uppercase;color:var(--mut)}
 .row{display:flex;align-items:center;gap:8px;padding:4px 0}.row label{display:flex;align-items:center;gap:7px;flex:1;cursor:pointer}
 .sw{width:12px;height:12px;border-radius:3px;display:inline-block;flex:none}
 .mutd{color:var(--mut);font-size:11.5px}.sp{height:10px}
 input[type=range]{width:130px}
 .chips{display:flex;flex-wrap:wrap;gap:5px}.chip{padding:2px 9px;border:1px solid var(--line);border-radius:999px;background:#fff;cursor:pointer;font-size:12px}
 .chip.on{background:var(--acc);border-color:var(--acc);color:#fff}.chip.pk{background:#f1e4fa;border-color:#b04fe0;color:#5a2078}
 .list{max-height:170px;overflow:auto;border:1px solid var(--line);border-radius:8px;padding:2px 8px}
 #wrap{flex:1;min-height:0;display:grid;position:relative}
 .pn{position:relative;min-height:0;border-left:1px solid var(--line)}.pn:first-child{border-left:0}
 .mv{position:absolute;inset:0}
 .ttl{position:absolute;left:12px;top:10px;z-index:5;color:#fff;text-shadow:0 1px 3px rgba(0,0,0,.6);pointer-events:none}
 .light .ttl{color:var(--ink);text-shadow:none}.ttl b{font-size:14px}.ttl span{display:block;font-size:11.5px;opacity:.85}
 .fbn{position:absolute;left:50%;bottom:14px;transform:translateX(-50%);z-index:6;padding:6px 14px;border-radius:999px;font-weight:700;font-size:12.5px;display:none;background:rgba(20,22,26,.86);color:#fff;border:2px solid;pointer-events:none;white-space:nowrap}
 #key{position:absolute;left:12px;bottom:12px;z-index:20;width:210px;padding:10px 12px;border-radius:10px;background:rgba(20,22,26,.78);color:#f1f1ef;font-size:12px}
 #key h5{margin:8px 0 5px;font-size:10px;letter-spacing:.08em;text-transform:uppercase;color:#aeb1ba}#key h5:first-child{margin-top:0}
 #key .it{display:flex;align-items:center;gap:7px;padding:2px 0;cursor:default;word-break:break-all}#key .it.tg{cursor:pointer}#key .it.off{opacity:.4;text-decoration:line-through}
 #key .sq{width:11px;height:11px;border-radius:3px;flex:none}
 .light #key{background:rgba(255,255,255,.92);color:var(--ink);box-shadow:0 1px 6px rgba(0,0,0,.15)}.light #key h5{color:var(--mut)}
 .msp-viewport-top-left-controls{display:none !important}
 #busy{position:absolute;right:56px;top:12px;z-index:25;font-size:12px;color:#fff;background:rgba(20,22,26,.7);padding:3px 10px;border-radius:999px;display:none}
</style>
<div id="root">
 <div id="bar">
  <div class="grp"><span class="gl">Surface</span>
   <div class="seg" id="sf"><button data-m="sas" class="on" style="--c:#2a78d6">SAS</button><button data-m="off" style="--c:#2a78d6">Off</button></div></div>
  <div class="grp">
   <button class="tool" data-p="pPk">Pockets ▾</button>
   <button class="tool" data-p="pRes">Residues ▾</button>
   <button class="tool" data-p="pView">View ▾</button>
   <button class="tool" id="bReset">Reset view</button>
   <button class="tool" id="bFs" title="Full-screen explorer (Esc to exit)" style="margin-left:auto;font-size:15px;padding:2px 9px">⛶</button></div>

  <div class="pop" id="pPk"><h4>Predicted pockets</h4><div id="pkrows"></div>
    <div class="sp"></div><button class="tool" id="bAll">Show all</button> <button class="tool" id="bNone">Hide all</button></div>

  <div class="pop" id="pRes" style="width:330px;max-width:330px"><h4>Highlight</h4>
    <div id="hl"><div class="row"><label><input type="checkbox" data-k="ligand">Ligand</label></div>
    <div class="row"><label><input type="checkbox" data-k="target">Target residue</label></div></div>
    <div class="sp"></div><h4>Residue colour</h4>
    <div class="seg" id="rc"><button data-v="pocket" class="on" style="--c:#2a78d6">By pocket</button><button data-v="type" style="--c:#2a78d6">By type</button><button data-v="score" style="--c:#2a78d6">By score</button></div>
    <div class="sp"></div><h4>Residues near the visible pockets</h4>
    <div class="row"><span class="mutd" style="width:60px">within</span><input id="sNbp" type="range" min="0" max="12" step="0.5" value="0"><b id="vNbp" style="width:46px">off</b></div>
    <div class="chips" id="types"></div><div class="sp"></div>
    <div class="list" id="nblist"></div>
    <div class="sp"></div><h4>Highlighted residues</h4><div class="chips" id="picks"></div><div class="sp"></div><button class="tool" id="bClr">Clear highlighted</button></div>

  <div class="pop" id="pView"><h4>Theme</h4>
    <div class="seg" id="th"><button data-t="dark" class="on" style="--c:#2a78d6">Dark</button><button data-t="light" style="--c:#2a78d6">Light</button></div>
    <div class="sp"></div><h4>Surface</h4>
    <div class="row"><span class="mutd" style="width:90px">Opacity</span><input id="sOp" type="range" min="15" max="100" value="100"></div>
    <div class="row"><span class="mutd" style="width:90px">Zoom speed</span><input id="sZoom" type="range" min="1" max="20" value="20"></div>
    <div class="sp"></div><h4>Show</h4>
    <div class="row"><label><input id="cEns" type="checkbox">All Lacuna conformations (ensemble)</label></div>
    <div class="row"><label><input id="cPol" type="checkbox">Polar contacts</label></div>
    <div class="row"><label><input id="cWh" type="checkbox" checked>Warhead markers</label></div>
    <div class="row"><label><input id="cRef" type="checkbox" checked>Experimental reference ligand</label></div>
    <div class="row"><label><input id="cLab" type="checkbox">Residue labels (pocket / highlighted)</label></div>
    <div class="row"><label><input id="cNbL" type="checkbox" checked>Residues around ligand</label></div>
    <div class="sp"></div><button class="tool" id="bPng">Save PNG</button>
    <div class="sp"></div><div class="mutd">Drag = rotate · right-drag = move · wheel = zoom · hover = name · Mol* tools: icons at top right of the view</div></div>
 </div>
 <div id="wrap"><div id="key"></div><div id="busy">updating…</div></div>
</div>
<script src="__CDN__"></script>
<script>
const SPEC = __SPEC__;
const wrap = document.getElementById('wrap'); wrap.style.gridTemplateColumns = 'repeat(' + SPEC.panels.length + ', 1fr)';
const S = {fade:1.0, hl:{ligand:false,pocket:false,target:false}, surf:'sas', ens:false, sop:1, pol:false, wh:true, ref:true, lab:false, nbl:true, rcol:'pocket', nbp:0, types:{}, dark:SPEC.dark, pkVis:{}, picked:{}};
const HLC = {ligand:'#00a9c0', pocket:'#e6781f', target:'#f2c200'}, HLN = {ligand:'Ligand', pocket:'Pocket', target:'Target residue'};
const PICK = '#b04fe0', TARGET = '#f2c200';
const SAS_DARK = '#e2e2ff', SAS_LIGHT = '#d8d8d2';
const CLASS_COL = {hydrophobic:'#e0a030', polar:'#2fb5a0', positive:'#3b6fd6', negative:'#d94a4a', 'Gly / Pro':'#9b8fd9', Cys:'#e87ba4'};
const TYPES = ['CYS','SER','LYS','TYR','HIS','THR','ASP','GLU','ARG'];
const POSEHEX = {cyan:'#00d5ea', magenta:'#ff4fd8', green:'#8fdc2b', orange:'#ff9d2e', purple:'#a47bff', white:'#e6e6e6'};
const P = []; window.__P = P; let camReset = true, queue = Promise.resolve(), timer = null, syncing = false;

// ---------- helpers
const resSel = spec => { const [ch,rest] = spec.split(':'); return {auth_asym_id: ch, auth_seq_id: parseInt(rest.replace(/[^0-9-]/g,''))}; };
const tint = (hex, t) => '#' + [1,3,5].map(i => { const v = parseInt(hex.substr(i,2),16); return Math.round(v + (255 - v) * t).toString(16).padStart(2,'0'); }).join('');   // mix toward white
const shade = (hex, t) => '#' + [1,3,5].map(i => Math.round(parseInt(hex.substr(i,2),16) * (1 - t)).toString(16).padStart(2,'0')).join('');   // mix toward black
const keyOf = spec => { const q = resSel(spec); return q.auth_asym_id+':'+q.auth_seq_id; };

// highlight "cloud": equal-sized translucent spheres on each highlighted residue's centroid; they overlap into one smooth, proportional cloud (camera is never moved)
function resCentroids(R){ if (R._rc) return R._rc; const acc = {}; R.atoms.forEach(a => { const k = a.chain+':'+a.resi; (acc[k] = acc[k] || [0,0,0,0]); acc[k][0]+=a.x; acc[k][1]+=a.y; acc[k][2]+=a.z; acc[k][3]++; });
  R._rc = {}; Object.keys(acc).forEach(k => { const v = acc[k]; R._rc[k] = [v[0]/v[3], v[1]/v[3], v[2]/v[3]]; }); return R._rc; }
function cloud(R, keys, color, radius){ const rc = resCentroids(R); return keys.filter(k => rc[k]).map(k => N('primitive', {kind:'ellipsoid', center: rc[k], radius, color})); }
const sel = k => { const [c,r] = k.split(':'); return {auth_asym_id:c, auth_seq_id:parseInt(r)}; };
const visible = pk => S.pkVis[pk.id] !== false;
const anyHL = () => S.hl.ligand || S.hl.pocket || S.hl.target;
// keep highlights / sliders / camera across Streamlit reruns (ticking a row, switching lane, resizing the window all rebuild this iframe)
const VKEY = 'pscreen_view2_' + (SPEC.panels[0].receptor.length) + '_' + (SPEC.panels[0].title || '');
let savedCam = null;
try { const sv = JSON.parse(sessionStorage.getItem(VKEY) || 'null'); if (sv) { if (sv.hl) { S.hl = {ligand: !!sv.hl.ligand, pocket: false, target: !!sv.hl.target}; } savedCam = sv.cam || null; } } catch(e) {}
const saveView = (cam) => { try { const o = JSON.parse(sessionStorage.getItem(VKEY) || '{}'); o.hl = S.hl; if (cam) o.cam = cam; sessionStorage.setItem(VKEY, JSON.stringify(o)); } catch(e) {} };
const visPoses = R => R.poses.filter(o => !o.p.site || S.pkVis[o.p.site] !== false);
const N = (kind, params, children) => { const n = {kind}; if (params !== undefined) n.params = params; if (children && children.length) n.children = children; return n; };
function resColor(pk, r){ if (S.rcol==='score' && pk.res_score && pk.res_score[r]) return pk.res_score[r]; if (S.rcol==='type') return (pk.res_class && pk.res_class[r]) || pk.color; return pk.color; }
function parsePdb(text){ const a = []; text.split('\n').forEach(l => { if (l.startsWith('ATOM') || l.startsWith('HETATM')) a.push({chain:l[21], resi:parseInt(l.slice(22,26)), resn:l.slice(17,20).trim(), x:+l.slice(30,38), y:+l.slice(38,46), z:+l.slice(46,54)}); }); return a; }
function parseSdf(text){ const L = text.split('\n'); const n = parseInt(L[3].slice(0,3)); const a = []; for (let i = 0; i < n; i++) { const l = L[4+i]; a.push({x:+l.slice(0,10), y:+l.slice(10,20), z:+l.slice(20,30)}); } return a; }
function centroids(atoms){ const acc = {}; atoms.forEach(a => { const k = a.chain+':'+a.resi, q = acc[k] = acc[k] || {x:0,y:0,z:0,n:0,resn:a.resn}; q.x+=a.x; q.y+=a.y; q.z+=a.z; q.n++; }); Object.keys(acc).forEach(k => { const q = acc[k]; q.x/=q.n; q.y/=q.n; q.z/=q.n; }); return acc; }
const blob = t => URL.createObjectURL(new Blob([t]));


// Mol*'s wheel zoom is proportional to the current camera distance, so once zoomed in close it crawls back out (and can end up inside the protein).
// Own zoom: step proportional to distance but with a floor, and the camera cannot come closer than MIN_D to its target.
const MIN_D = 9;
// Whatever moved the camera (wheel, pinch, pan, saved view), it is pushed back out if it ended up inside the protein body (flat grey void)
const CLEAR2 = 72;      // (8.5 Å)^2 from any atom centre ≈ 5 Å clear of the surface
// Mol*'s default clipping (radius 48 %) slices the front of the surface and fog dims it: both show up as grey holes/regions while rotating. Mol* (or an MVS reload) resets these, so they are re-asserted whenever the camera moves.
function ensureDepth(R){
  const c = R.viewer.plugin.canvas3d; if (!c) return; const p = c.props;
  if (p.cameraClipping.radius !== 0 || p.cameraClipping.far || p.cameraFog.name !== 'off') c.setProps({cameraFog:{name:'off', params:{}}, cameraClipping:{radius:0, far:false, minNear:1}});
}
function keepOutside(R){
  ensureDepth(R);
  const cam = R.viewer.plugin.canvas3d && R.viewer.plugin.canvas3d.camera; if (!cam || R._ko) return;
  const s = cam.getSnapshot(); let p = s.position.slice(), moved = false;
  for (let it = 0; it < 8; it++) {
    let best = null, bd = CLEAR2; for (const a of R.atoms) { const d = (a.x-p[0])**2 + (a.y-p[1])**2 + (a.z-p[2])**2; if (d < bd) { bd = d; best = a; } }
    if (!best) break; const v = [p[0]-best.x, p[1]-best.y, p[2]-best.z], l = Math.hypot(...v) || 1, k = Math.sqrt(CLEAR2) * 1.02 / l; p = [best.x+v[0]*k, best.y+v[1]*k, best.z+v[2]*k]; moved = true;
  }
  if (moved) { R._ko = true; try { cam.setState({position: p}, 0); } finally { setTimeout(() => { R._ko = false; }, 0); } }
}
function installWheelZoom(host, viewer){
  host.addEventListener('wheel', ev => {
    const cam = viewer.plugin.canvas3d && viewer.plugin.canvas3d.camera; if (!cam) return;
    ev.preventDefault(); ev.stopPropagation(); ev.stopImmediatePropagation();
    const s = cam.getSnapshot(), t = s.target, v = [0,1,2].map(i => s.position[i] - t[i]), d = Math.hypot(...v); if (d < 1e-6) return;
    const spd = parseFloat(document.getElementById('sZoom').value) / 20;                       // slider 1-20, default 20 = full speed
    const dy = Math.max(-100, Math.min(100, ev.deltaMode === 1 ? ev.deltaY * 33 : ev.deltaY));
    let nd = d * Math.exp(dy * 0.0016 * spd); if (dy > 0) nd = Math.max(nd, d + 2.5 * spd * dy / 100);   // zooming out always moves at least a few Å per tick
    nd = Math.max(MIN_D, nd);
    // never let the camera enter the protein body: from inside, the back faces of the surface show as a flat grey void
    const R = P.find(r => r.viewer === viewer);
    if (R && nd < d) { const inside = f => { const q = [0,1,2].map(i => t[i] + v[i] * f / d); return R.atoms.some(a => (a.x-q[0])**2 + (a.y-q[1])**2 + (a.z-q[2])**2 < 72); };   // 8.5 Å from any atom ≈ 5 Å clear of the surface, so the near plane never slices it
      if (inside(nd)) { let lo = nd, hi = d; for (let i = 0; i < 12; i++) { const m = (lo + hi) / 2; if (inside(m)) lo = m; else hi = m; } nd = hi; } }
    const k = nd / d;
    cam.setState({position:[0,1,2].map(i => t[i] + v[i] * k)}, 0);   // radius (clipping range) must NOT shrink with the zoom, or the near plane slices the surface into grey holes
  }, {capture:true, passive:false});
}
// ---------- panels
(async () => {
  for (const [k, pd] of SPEC.panels.entries()) {
    const el = document.createElement('div'); el.className = 'pn'; el.innerHTML = '<div class="mv" id="mv'+k+'"></div><div class="ttl"><b>'+pd.title+'</b><span>'+pd.subtitle+'</span></div><div class="fbn"></div>'; wrap.appendChild(el);
    const viewer = await molstar.Viewer.create('mv'+k, {layoutIsExpanded:false, layoutShowControls:false, layoutShowSequence:false, layoutShowLog:false, layoutShowLeftPanel:false, viewportShowExpand:false, viewportShowSelectionMode:false, viewportShowAnimation:false});
    installWheelZoom(document.getElementById('mv'+k), viewer);
    const atoms = parsePdb(pd.receptor);
    P.push({k, el, pd, viewer, atoms, cent: centroids(atoms), urlRec: blob(pd.receptor), urlRef: pd.reference ? blob(pd.reference.data) : null, urlEns: pd.ensemble ? blob(pd.ensemble) : null,
            poses: pd.poses.map(p => ({p, url: blob(p.data), atoms: parseSdf(p.data)})), _nbk:null, _nb:[]});
  }
  document.getElementById('wrap').appendChild(document.getElementById('key')); document.getElementById('wrap').appendChild(document.getElementById('busy'));
  P[0].pd.pockets.forEach(pk => S.pkVis[pk.id] = true);
  (SPEC.panels[0].picked || []).forEach(k => { const [c,r] = k.split(':'); const m = r.match(/([A-Z]{3})?(-?\d+)/); S.picked[c+':'+m[2]] = m[1]||''; });
  buildUI(); schedule(!savedCam);
  P.forEach(R => { try { R.viewer.plugin.canvas3d.camera.stateChanged.subscribe(() => keepOutside(R)); } catch(e) {} });
  try { let tm = null; P[0].viewer.plugin.canvas3d.camera.stateChanged.subscribe(() => { clearTimeout(tm); tm = setTimeout(() => { try { saveView(P[0].viewer.plugin.canvas3d.camera.getSnapshot()); } catch(e) {} }, 400); }); } catch(e) {}
  if (P.length > 1) { try { P.forEach((a,i) => a.viewer.plugin.canvas3d.camera.stateChanged.subscribe(() => { if (syncing) return; syncing = true; const snap = a.viewer.plugin.canvas3d.camera.getSnapshot(); P.forEach((b,j) => { if (j!==i) b.viewer.plugin.canvas3d.camera.setState(snap, 0); }); syncing = false; })); } catch(e) {} }
})();

// ---------- the scene, declaratively (MolViewSpec)
function nearbyPocketResidues(R){
  if (!(S.nbp > 0)) return [];
  const pk = visPocketKeys(R), tf = Object.keys(S.types).filter(t => S.types[t]).join(','), key = [...pk].join('|')+'/'+S.nbp+'/'+tf; if (R._nbk === key) return R._nb;
  const pkAtoms = R.atoms.filter(a => pk.has(a.chain+':'+a.resi)), best = new Map(), resn = {};
  R.atoms.forEach(a => { const k = a.chain+':'+a.resi; if (pk.has(k)) return; let d2 = Infinity; for (const b of pkAtoms) { const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z, q = dx*dx+dy*dy+dz*dz; if (q < d2) d2 = q; }
    if (!best.has(k) || d2 < best.get(k)) { best.set(k, d2); resn[k] = a.resn; } });
  let out = [...best.entries()].filter(([k,d2]) => d2 <= S.nbp*S.nbp).map(([k,d2]) => ({key:k, resn:resn[k], d:Math.sqrt(d2)}));
  if (tf) out = out.filter(o => S.types[o.resn]); out.sort((a,b) => a.d - b.d); R._nbk = key; R._nb = out; return out;
}
function ownedResidues(R){ const seen = {}, out = {}, ovl = new Set(), vis = R.pd.pockets.filter(visible);
  vis.forEach(pk => pk.residues.forEach(r => { const k = keyOf(r); (seen[k] = seen[k] || new Set()).add(pk.model || pk.id); }));
  Object.keys(seen).forEach(k => { if (seen[k].size > 1) ovl.add(k); });
  vis.forEach(pk => { out[pk.id] = pk.residues.filter(r => !ovl.has(keyOf(r))); }); out.__overlap = [...ovl]; return out; }
function visPocketKeys(R){ const s = new Set(); R.pd.pockets.forEach(pk => { if (visible(pk)) pk.residues.forEach(r => s.add(keyOf(r))); }); return s; }
function aroundLigand(R){ const vp = visPoses(R), out = new Set(); if (!vp.length) return out; vp.forEach(o => o.atoms.forEach(la => { R.atoms.forEach(a => { const dx=a.x-la.x, dy=a.y-la.y, dz=a.z-la.z; if (dx*dx+dy*dy+dz*dz <= 16) out.add(a.chain+':'+a.resi); }); })); return out; }
const ELEM = [['N','#6f86ff'],['O','#ff5a5a'],['S','#f0d24a'],['H','#f0f0f0']];
const sticks = (carbon, size) => N('representation', {type:'ball_and_stick', size_factor:size}, [N('color',{color:carbon}), ...ELEM.map(([e,c]) => N('color',{color:c, selector:{type_symbol:e}}))]);

function buildMVS(R){
  const hl = anyHL(), base = S.dark ? SAS_DARK : SAS_LIGHT, vp = visPoses(R), kids = [];
  kids.push(N('canvas', {background_color: S.dark ? '#14143a' : '#ffffff'}));
  // ----- receptor
  const rec = [], visKeys = visPocketKeys(R);
  const polyReps = [];
  if (S.surf === 'sas' && !S.ens) {
    const patches = [];
    const own = ownedResidues(R); R.pd.pockets.forEach(pk => { if (visible(pk) && own[pk.id].length) patches.push(N('color', {color: pk.color, selector: own[pk.id].map(resSel)})); });
    if (own.__overlap.length) patches.push(N('color', {color:'#c23fe0', selector: own.__overlap.map(k => { const [c,n] = k.split(':'); return {auth_asym_id:c, auth_seq_id:parseInt(n)}; })}));
    if (R.pd.target) patches.push(N('color', {color:TARGET, selector: resSel(R.pd.target)}));
    const pk2 = Object.keys(S.picked); if (pk2.length) patches.push(N('color', {color:PICK, selector: pk2.map(sel)}));
    // pocket highlight: single opaque surface, pocket in its own colour, background untouched (no overlay layers, no tint, no fade)
    polyReps.push(N('representation', {type:'surface'}, [N('color',{color:base}), N('opacity',{opacity: S.sop}), ...patches]));
  }
  if (polyReps.length) rec.push(N('component', {selector:'polymer'}, polyReps));
  // pocket residues (sticks; solid surface + labels while the pocket is highlighted)
  const own2 = ownedResidues(R);
  R.pd.pockets.forEach(pk => { if (!visible(pk) || !own2[pk.id].length) return;
    const reps = [sticks(pk.color, S.hl.pocket ? 0.7 : 0.4)];
    rec.push(N('component', {selector: own2[pk.id].map(resSel)}, reps));
    if (S.hl.pocket && S.lab) own2[pk.id].forEach(r => rec.push(N('component', {selector:resSel(r)}, [N('label',{text:r.split(':')[1]})])));
  });
  if (own2.__overlap.length) rec.push(N('component', {selector: own2.__overlap.map(k => { const [c,n] = k.split(':'); return {auth_asym_id:c, auth_seq_id:parseInt(n)}; })}, [sticks('#c23fe0', S.hl.pocket ? 0.7 : 0.45)]));
  if (S.nbl) { const ar = [...aroundLigand(R)].filter(k => !visKeys.has(k) && !S.picked[k]); if (ar.length) rec.push(N('component', {selector: ar.map(sel)}, [sticks(S.dark?'#d8d8e8':'#777', 0.3)])); }
  const nb = nearbyPocketResidues(R).filter(o => !S.picked[o.key]); if (nb.length) rec.push(N('component', {selector: nb.map(o => sel(o.key))}, [sticks(S.dark?'#d8d8e8':'#777', 0.3)]));
  Object.keys(S.picked).forEach(k => { const reps = [sticks(PICK, 0.7)]; rec.push(N('component', {selector: sel(k)}, reps)); if (S.lab) rec.push(N('component', {selector: sel(k)}, [N('label',{text:(S.picked[k]||'')+k.split(':')[1]})])); });
  let TGT_LABEL = null;
  if (R.pd.target) { const t = resSel(R.pd.target), reps = [sticks(TARGET, S.hl.target ? 1.1 : 0.6)];
    rec.push(N('component', {selector:t}, reps)); const ta = R.atoms.filter(a => a.chain === t.auth_asym_id && a.resi === t.auth_seq_id); if (ta.length) { const tc = [0,1,2].map(i => ta.reduce((m,a) => m + [a.x,a.y,a.z][i], 0) / ta.length); TGT_LABEL = {kind:'label', position:tc, text:R.pd.target.split(':')[1], label_size:3, label_color:'#ffd21f', label_offset:8}; } }
  // camera focus target
  // highlights never move the camera: the view stays where the user put it (only Reset / first load frames the scene)
  kids.push(N('download', {url:R.urlRec}, [N('parse', {format:'pdb'}, [N('structure', {type:'model'}, rec)])]));
  if (S.ens && R.urlEns && R.pd.ens_info) { const ei = R.pd.ens_info, ms = [];
    for (let i = 0; i < ei.n; i++) { const opens = ei.open.indexOf(i) >= 0, inp = i === 0, c = inp ? '#8d949e' : (opens ? (i === ei.docked ? '#f28c28' : '#3d7ff0') : '#9db6d6');
      ms.push(N('structure', {type:'model', model_index:i}, [N('component', {selector:'polymer'}, [N('representation', {type:'cartoon', size_factor: (opens || inp) ? 0.55 : 0.18}, [N('color',{color:c}), N('opacity',{opacity: (opens || inp) ? 0.95 : 0.4})])])])); }
    kids.push(N('download', {url:R.urlEns}, [N('parse', {format:'pdb'}, ms)])); }

  // ----- poses
  vp.forEach((o, i) => { const reps = [sticks(POSEHEX[o.p.color], S.hl.ligand ? 0.9 : 0.6)];
    if (S.hl.ligand) reps.push(N('representation', {type:'surface'}, [N('color',{color:'#00d5ea'}), N('opacity',{opacity:0.5})]));
    const comp = [...reps, N('tooltip', {text:o.p.label})];
    kids.push(N('download', {url:o.url}, [N('parse', {format:'sdf'}, [N('structure', {type:'model'}, [N('component', {selector:'all'}, comp)])])])); });
  if (R.urlRef && S.ref) kids.push(N('download', {url:R.urlRef}, [N('parse', {format:'pdb'}, [N('structure', {type:'model'}, [N('component', {selector:'all'}, [sticks('#5be37d', 0.4)])])])]));
  // ----- primitives: target distance, warhead markers, polar contacts
  const prim = [];
  vp.forEach(o => { const q = o.p;
    if (S.wh) { q.warhead_xyz.forEach(c => prim.push(N('primitive', {kind:'ellipsoid', center:c, radius:0.28, color:'#ffb066'}))); q.reactive_xyz.forEach(c => prim.push(N('primitive', {kind:'ellipsoid', center:c, radius:0.42, color:'#ff3b3b'}))); }
    if (S.pol) q.polar.forEach(c => prim.push(N('primitive', {kind:'distance_measurement', start:[c[0],c[1],c[2]], end:[c[3],c[4],c[5]], label_template:c[6].toFixed(1), color:'#ffe14d', label_color:'#ffe14d', radius:0.07, label_size:1.3})));
  });
  if (TGT_LABEL) prim.push(N('primitive', TGT_LABEL));
  const cl = [];
  if (S.hl.target && R.pd.target) cl.push(...cloud(R, [keyOf(R.pd.target)], TARGET, 4.6));
  if (cl.length) kids.push(N('primitives', {opacity: 0.5}, cl));
  if (prim.length) kids.push(N('primitives', {}, prim));
  // ----- camera (initial / reset): frame poses + pockets
  if (camReset) { const pts = []; vp.forEach(o => o.atoms.forEach(a => pts.push([a.x,a.y,a.z]))); visKeys.forEach(k => { const c = R.cent[k]; if (c) pts.push([c.x,c.y,c.z]); });
    if (pts.length) { const m = [0,1,2].map(i => pts.reduce((s,p) => s+p[i], 0)/pts.length), rad = Math.max(...pts.map(p => Math.hypot(p[0]-m[0], p[1]-m[1], p[2]-m[2])));
      const pc = [0,1,2].map(i => R.atoms.reduce((s,a) => s + [a.x,a.y,a.z][i], 0)/R.atoms.length);              // protein centre
      let d = [m[0]-pc[0], m[1]-pc[1], m[2]-pc[2]]; const L = Math.hypot(...d); d = L < 1 ? [0,0,1] : d.map(x => x/L);          // look at the pocket from outside the protein
      const dist = Math.max(34, rad*2.8), up = Math.abs(d[1]) > 0.95 ? [1,0,0] : [0,1,0];
      kids.push(N('camera', {target:m, position:[m[0]+d[0]*dist, m[1]+d[1]*dist, m[2]+d[2]*dist], up})); } }
  return {metadata:{version:'1', timestamp:new Date().toISOString()}, root:N('root', undefined, kids)};
}

async function renderAll(){
  document.getElementById('busy').style.display = 'block';
  for (const R of P) { try { await R.viewer.loadMvsData(JSON.stringify(buildMVS(R)), 'mvsj', {replaceExisting:true, keepCamera: !camReset}); } catch(e) { console.error('MVS error', e); } }
  camReset = false; document.getElementById('busy').style.display = 'none'; applyZoom(); P.forEach(ensureDepth);
  if (savedCam) { try { P.forEach(R => R.viewer.plugin.canvas3d.camera.setState(savedCam, 0)); } catch(e) {} savedCam = null; }
}
// Mol*'s hover box lists implementation details (blob URL, model, instance, entity); keep only residue / compound names
new MutationObserver(() => document.querySelectorAll('.msp-highlight-simple-row').forEach(r => { const t = r.textContent.trim(); r.style.display = (/^(blob:|data:|Model \d|Instance |Polymer \d|Unknown Entity)/.test(t)) ? 'none' : ''; })).observe(document.body, {childList:true, subtree:true, characterData:true});
function schedule(reset){ saveView(reset ? false : null); if (reset) { try { const o = JSON.parse(sessionStorage.getItem(VKEY) || '{}'); delete o.cam; sessionStorage.setItem(VKEY, JSON.stringify(o)); } catch(e) {} savedCam = null; camReset = true; } clearTimeout(timer); timer = setTimeout(() => { queue = queue.then(renderAll); }, 140); updateUI(); }
function applyZoom(){ try { const z = parseFloat(document.getElementById('sZoom').value) * 0.03; P.forEach(R => R.viewer.plugin.canvas3d.setProps({trackball:{zoomSpeed: z}, cameraFog:{name:'off', params:{}}, cameraClipping:{radius:0, far:false, minNear:1}})); } catch(e) {} }

// ---------- UI
function updateUI(){
  document.querySelectorAll('#hl input').forEach(b => { b.checked = !!S.hl[b.dataset.k]; });
  const act = Object.keys(S.hl).filter(k => S.hl[k]);
  P.forEach(R => { const b = R.el.querySelector('.fbn'); if (!act.length) { b.style.display='none'; return; } b.style.display='block'; b.style.borderColor = HLC[act[0]]; b.innerHTML = act.map(k => '<span style="color:'+HLC[k]+'">●</span> '+HLN[k]).join(' &nbsp;+&nbsp; '); });
  document.getElementById('root').classList.toggle('light', !S.dark); P.forEach(R => R.el.classList.toggle('light', !S.dark));
  if (!P.length) return; const R0 = P[0], vp = visPoses(R0); let h = '';
  const pks = R0.pd.pockets; if (pks.length) h += '<h5>Pockets</h5>' + pks.map(pk => '<div class="it tg '+(visible(pk)?'':'off')+'" data-pk="'+pk.id+'"><span class="sq" style="background:'+pk.color+'"></span>'+pk.id+'</div>').join('');
  const names = []; vp.forEach(o => { if (!names.find(n => n.l === o.p.label)) names.push({l:o.p.label, c:o.p.color, d:o.p.distance ? o.p.distance.text : ''}); });
  if (names.length) h += '<h5>Compounds</h5>' + names.map(n => '<div class="it"><span class="sq" style="background:'+POSEHEX[n.c]+'"></span>'+n.l+(n.d && R0.pd.target ? ' <span class="mutd" style="color:#ffd21f">· '+n.d+' to '+R0.pd.target.split(':')[1]+'</span>' : '')+'</div>').join('');
  const mk = [];
  if (R0.pd.target) mk.push([TARGET, 'Target '+R0.pd.target.split(':')[1]]);
  if (Object.keys(S.picked).length) mk.push([PICK, 'Highlighted residues']);
  if (R0.urlRef && S.ref) mk.push(['#5be37d', 'Experimental ligand']);
  if (ownedResidues(R0).__overlap.length) mk.push(['#c23fe0', 'Overlap (2+ approaches)']);
  if (S.wh && vp.some(o => o.p.reactive_xyz.length)) mk.push(['#ff3b3b', 'Reactive atom']);
  if (mk.length) h += '<h5>Markers</h5>' + mk.map(m => '<div class="it"><span class="sq" style="background:'+m[0]+'"></span>'+m[1]+'</div>').join('');
  if (S.rcol === 'type' && pks.some(visible)) h += '<h5>Residue type</h5>' + Object.keys(CLASS_COL).map(k => '<div class="it"><span class="sq" style="background:'+CLASS_COL[k]+'"></span>'+k+'</div>').join('');
  if (S.ens && R0.pd.ens_info) { const ei = R0.pd.ens_info; h += '<h5>Lacuna conformations ('+ei.n+')</h5>' + [['#8d949e','Input structure (0)'], ['#f28c28','Docked: conformer '+ei.docked], ['#3d7ff0','Other pocket-opening conformers'], ['#9db6d6','No pocket detected']].map(m => '<div class="it"><span class="sq" style="background:'+m[0]+'"></span>'+m[1]+'</div>').join(''); }
  const kk = document.getElementById('key'); kk.innerHTML = h; kk.style.display = h ? 'block' : 'none';
  kk.querySelectorAll('[data-pk]').forEach(e => e.onclick = () => { S.pkVis[e.dataset.pk] = !visible({id:e.dataset.pk}); syncPk(); schedule(false); });
}
function syncPk(){ document.querySelectorAll('#pkrows input').forEach(i => i.checked = visible({id:i.dataset.id})); }
function refreshPicks(){
  const pk = document.getElementById('picks'); pk.innerHTML = Object.keys(S.picked).length ? Object.keys(S.picked).map(k => '<span class="chip pk">'+(S.picked[k]||'')+k.split(':')[1]+' ('+k.split(':')[0]+') <b data-k="'+k+'" style="cursor:pointer">×</b></span>').join('') : '<span class="mutd">none</span>';
  pk.querySelectorAll('b').forEach(b => b.onclick = () => { delete S.picked[b.dataset.k]; refreshPicks(); schedule(false); });
  const list = document.getElementById('nblist'); if (!P.length) return; const nb = nearbyPocketResidues(P[0]);
  if (!(S.nbp > 0)) { list.innerHTML = '<div class="mutd" style="padding:8px 0">Increase “within” to list residues near the visible pockets. Tick a residue to highlight it.</div>'; return; }
  list.innerHTML = nb.length ? nb.slice(0, 80).map(o => '<div class="row"><label><input type="checkbox" data-k="'+o.key+'" '+(S.picked[o.key]!==undefined?'checked':'')+'>'+o.resn+o.key.split(':')[1]+' <span class="mutd">chain '+o.key.split(':')[0]+' · '+o.d.toFixed(1)+' Å</span></label></div>').join('') : '<div class="mutd" style="padding:8px 0">No residues of the chosen types within this distance.</div>';
  list.querySelectorAll('input').forEach(i => i.onchange = () => { const o = nb.find(x => x.key === i.dataset.k); if (i.checked) S.picked[i.dataset.k] = o.resn; else delete S.picked[i.dataset.k]; refreshPicks(); schedule(false); });
}
function buildUI(){
  document.querySelectorAll('#hl input').forEach(b => b.onchange = () => { S.hl[b.dataset.k] = b.checked; schedule(false); });
  document.querySelectorAll('#sf button').forEach(b => b.onclick = () => { S.surf = b.dataset.m; document.querySelectorAll('#sf button').forEach(x => x.classList.toggle('on', x===b)); schedule(false); });
  document.getElementById('bReset').onclick = () => { S.hl = {ligand:false,pocket:false,target:false}; schedule(true); };
  document.querySelectorAll('.tool[data-p]').forEach(b => b.onclick = (e) => { e.stopPropagation(); const p = document.getElementById(b.dataset.p), open = p.classList.contains('open');
    document.querySelectorAll('.pop').forEach(x => x.classList.remove('open')); document.querySelectorAll('.tool[data-p]').forEach(x => x.classList.remove('open'));
    if (!open) { p.classList.add('open'); b.classList.add('open'); const r = b.getBoundingClientRect(); p.style.left = Math.max(8, Math.min(r.left, window.innerWidth - p.offsetWidth - 10)) + 'px'; } });
  document.addEventListener('click', (e) => { if (!e.target.closest('.pop')) { document.querySelectorAll('.pop').forEach(x => x.classList.remove('open')); document.querySelectorAll('.tool[data-p]').forEach(x => x.classList.remove('open')); } });
  const pkr = document.getElementById('pkrows');
  P[0].pd.pockets.forEach(pk => { const d = document.createElement('div'); d.className = 'row';
    d.innerHTML = '<label><input type="checkbox" checked data-id="'+pk.id+'"><span class="sw" style="background:'+pk.color+'"></span><b>'+pk.id+'</b> <span class="mutd">'+pk.residues.length+' residues'+(pk.source ? ' · '+pk.source.replace(/_/g,' ') : '')+'</span></label>'; pkr.appendChild(d);
    d.querySelector('input').onchange = (e) => { S.pkVis[pk.id] = e.target.checked; refreshPicks(); schedule(false); }; });
  document.getElementById('bAll').onclick = () => { P[0].pd.pockets.forEach(p => S.pkVis[p.id] = true); syncPk(); refreshPicks(); schedule(false); };
  document.getElementById('bNone').onclick = () => { P[0].pd.pockets.forEach(p => S.pkVis[p.id] = false); syncPk(); refreshPicks(); schedule(false); };
  document.querySelectorAll('#rc button').forEach(b => b.onclick = () => { S.rcol = b.dataset.v; document.querySelectorAll('#rc button').forEach(x => x.classList.toggle('on', x===b)); schedule(false); });
  const tw = document.getElementById('types'); TYPES.forEach(t => { const c = document.createElement('span'); c.className = 'chip'; c.textContent = t; c.onclick = () => { S.types[t] = !S.types[t]; c.classList.toggle('on', S.types[t]); P.forEach(R => R._nbk = null); refreshPicks(); schedule(false); }; tw.appendChild(c); });
  document.getElementById('sNbp').oninput = (e) => { S.nbp = parseFloat(e.target.value); document.getElementById('vNbp').textContent = S.nbp > 0 ? S.nbp.toFixed(1)+' Å' : 'off'; refreshPicks(); schedule(false); };
  document.getElementById('bClr').onclick = () => { S.picked = {}; refreshPicks(); schedule(false); };
  document.querySelectorAll('#th button').forEach(b => b.onclick = () => { S.dark = b.dataset.t === 'dark'; document.querySelectorAll('#th button').forEach(x => x.classList.toggle('on', x===b)); schedule(false); });
  const bindc = (id, key) => { document.getElementById(id).onchange = (e) => { S[key] = e.target.checked; schedule(false); }; };
  bindc('cEns','ens'); bindc('cPol','pol'); bindc('cWh','wh'); bindc('cRef','ref'); bindc('cLab','lab'); bindc('cNbL','nbl');
  document.getElementById('sOp').oninput = (e) => { S.sop = parseFloat(e.target.value)/100; schedule(false); };
  document.getElementById('sZoom').oninput = applyZoom;
  document.getElementById('bPng').onclick = () => { try { P[0].viewer.plugin.helpers.viewportScreenshot.download(); } catch(e) { alert('Use the camera icon at the top right of the view.'); } };
  refreshPicks();
}

(function(){ const f = window.frameElement, P = window.parent; if (!f || !P) return; let full = false, base = null;
  const fit = () => { if (full) return; const m = P.document.querySelector('[data-testid="stMain"]') || P.document.documentElement;
    const top = f.getBoundingClientRect().top - m.getBoundingClientRect().top + m.scrollTop; f.style.height = Math.max(360, P.innerHeight - top - 32) + 'px'; };
  fit(); P.addEventListener('resize', fit); setTimeout(fit, 600); setTimeout(fit, 1800);
  const set = (on) => { full = on; const b = document.getElementById('bFs');
    if (on) { base = f.getAttribute('style') || ''; f.style.cssText = base + ';position:fixed;inset:0;width:100vw;height:100vh;z-index:999999;background:#fff;border:0'; b.textContent = '✕ Exit'; }
    else { f.style.cssText = base || ''; b.textContent = '⛶'; fit(); }
    setTimeout(() => window.dispatchEvent(new Event('resize')), 60); };
  document.getElementById('bFs').onclick = () => set(!full);
  document.addEventListener('keydown', e => { if (e.key === 'Escape' && full) set(false); });
})();
</script>
"""
