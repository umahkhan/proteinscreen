"""Export what the 3D explorer shows as a PyMOL script bundle (.pml + the structure files it loads), in the explorer's visual design.

Input is the explorer panel dict (pscreen.ui.viewer3d.panel), so the export always matches the current view: same pocket colours,
compound colours, element colours, target-residue highlight, warhead / reactive-atom markers and the warhead→target distance.
Pure text generation; PyMOL is only needed to open the result (cd into the unzipped folder, then `@pscreen_scene.pml`).
"""
from __future__ import annotations
import io, re, zipfile

BASE_DARK, BASE_LIGHT = "#e2e2ff", "#d8d8d2"
BG_DARK, BG_LIGHT = "#14143a", "#ffffff"
TARGET, OVERLAP = "#f2c200", "#c23fe0"
ELEM = {"N": "#6f86ff", "O": "#ff5a5a", "S": "#f0d24a"}
POSE = {"cyan": "#00d5ea", "magenta": "#ff4fd8", "green": "#8fdc2b", "orange": "#ff9d2e", "purple": "#a47bff", "white": "#e6e6e6"}
WARHEAD, REACTIVE, DIST = "#ffb066", "#ff3b3b", "#ffe14d"


def _rgb(h: str) -> str:
    h = h.lstrip("#"); return "[" + ", ".join(f"{int(h[i:i + 2], 16) / 255:.3f}" for i in (0, 2, 4)) + "]"


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(s)).strip("_") or "x"


def _sel(residues) -> str:
    """'A:LEU601' / 'A:601' ... -> PyMOL selection."""
    parts = []
    for r in residues:
        m = re.match(r"^(\w+):(?:[A-Za-z]{3})?(-?\d+)", str(r))
        if m: parts.append(f"(chain {m.group(1)} and resi {m.group(2)})")
    return " or ".join(parts)


def build(panel: dict, dark: bool = True) -> dict[str, str]:
    """Returns {filename: text}: pscreen_scene.pml, receptor.pdb, ligands/*.mol, [reference_ligand.pdb], README.txt."""
    files: dict[str, str] = {"receptor.pdb": panel["receptor"]}
    L: list[str] = []
    add = L.append
    colors: dict[str, str] = {}

    def col(hexc: str) -> str:
        n = "c_" + hexc.lstrip("#").lower(); colors[n] = hexc; return n

    base = BASE_DARK if dark else BASE_LIGHT
    cb, ct, co = col(base), col(TARGET), col(OVERLAP)
    pockets = [p for p in panel.get("pockets", []) if p.get("residues")]
    # residues claimed by more than one pocket get the overlap colour (as in the explorer)
    owners: dict[str, int] = {}
    for p in pockets:
        for r in set(p["residues"]): owners[r] = owners.get(r, 0) + 1
    overlap = [r for r, n in owners.items() if n > 1]

    add(f"# pscreen → PyMOL · {panel.get('title', '')} · {panel.get('subtitle', '')}")
    add("# Run from this folder: cd into it in PyMOL, then  @pscreen_scene.pml   (scenes: F1 overview · F2 pockets · F3 compounds · F4 covalent site)")
    add("reinitialize"); add("load receptor.pdb, receptor")
    body: list[str] = []
    b = body.append
    b("bg_color c_bg")
    b("hide everything"); b("set ray_opaque_background, on"); b("set antialias, 2"); b("set surface_quality, 1"); b("set cartoon_fancy_helices, 1")
    b("set stick_radius, 0.16"); b("set sphere_scale, 1.0"); b("set label_font_id, 7"); b("set label_size, 14"); b("set label_color, white" if dark else "set label_color, black")
    b("set transparency_mode, 2")

    # protein: base surface object (solid in overview, faded in the pocket scene)
    b("create prot, receptor and polymer"); b("show surface, prot"); b(f"color {cb}, prot"); b("set surface_color, " + cb + ", prot")
    b("delete receptor")

    pk_objs = []
    for i, p in enumerate(pockets):
        name = "pocket_" + _safe(p["id"]); pk_objs.append(name); cc = col(p.get("color") or "#e6935c")
        own = [r for r in p["residues"] if r not in overlap]
        if own:
            b(f"create {name}, prot and ({_sel(own)})")
            b(f"show surface, {name}"); b(f"show sticks, {name}"); b(f"color {cc}, {name}"); b(f"set surface_color, {cc}, {name}")
            for e, h in ELEM.items(): b(f"color {col(h)}, {name} and elem {e}")
        # pocket centre label
        if p.get("center"):
            cx, cy, cz = p["center"]; b(f"pseudoatom lbl_{name}, pos=[{cx:.3f}, {cy:.3f}, {cz:.3f}], label={p['id']}"); b(f"hide everything, lbl_{name}")
            b(f"label lbl_{name}, '{p['id']}'")
    if overlap:
        b(f"create pocket_overlap, prot and ({_sel(overlap)})"); pk_objs.append("pocket_overlap")
        b("show surface, pocket_overlap"); b("show sticks, pocket_overlap"); b(f"color {co}, pocket_overlap"); b(f"set surface_color, {co}, pocket_overlap")
        for e, h in ELEM.items(): b(f"color {col(h)}, pocket_overlap and elem {e}")

    # covalent / target residue of interest
    tgt = panel.get("target")
    if tgt:
        ts = _sel([tgt]); lab = tgt.split(":")[1]
        b(f"create target_res, prot and ({ts})"); b("show sticks, target_res"); b("set stick_radius, 0.24, target_res"); b(f"color {ct}, target_res")
        for e, h in ELEM.items(): b(f"color {col(h)}, target_res and elem {e}")
        b(f"label target_res and name CA, '{lab}'"); b("set label_size, 18, target_res")

    # compounds
    cmp_objs, mark_objs, dist_objs, polar_objs = [], [], [], []
    for i, p in enumerate(panel.get("poses", [])):
        nm = "cmpd_" + _safe(p["label"]); files[f"ligands/{nm}.mol"] = p["data"]; cmp_objs.append(nm)
        cc = col(POSE.get(p.get("color"), "#00d5ea"))
        b(f"load ligands/{nm}.mol, {nm}"); b(f"show sticks, {nm}"); b(f"set stick_radius, 0.22, {nm}"); b(f"color {cc}, {nm} and elem C")
        for e, h in ELEM.items(): b(f"color {col(h)}, {nm} and elem {e}")
        for k, (x, y, z) in enumerate(p.get("warhead_xyz") or []):
            mk = f"warhead_{nm}_{k}"; b(f"pseudoatom {mk}, pos=[{x:.3f}, {y:.3f}, {z:.3f}]"); b(f"show spheres, {mk}"); b(f"set sphere_scale, 0.28, {mk}"); b(f"color {col(WARHEAD)}, {mk}"); mark_objs.append(mk)
        for k, (x, y, z) in enumerate(p.get("reactive_xyz") or []):
            mk = f"reactive_{nm}_{k}"; b(f"pseudoatom {mk}, pos=[{x:.3f}, {y:.3f}, {z:.3f}]"); b(f"show spheres, {mk}"); b(f"set sphere_scale, 0.42, {mk}"); b(f"color {col(REACTIVE)}, {mk}"); mark_objs.append(mk)
        d = p.get("distance")
        if d:
            a, c2 = d["a"], d["b"]
            b(f"pseudoatom da_{nm}, pos=[{a[0]:.3f}, {a[1]:.3f}, {a[2]:.3f}]"); b(f"pseudoatom db_{nm}, pos=[{c2[0]:.3f}, {c2[1]:.3f}, {c2[2]:.3f}]")
            b(f"distance dist_{nm}, da_{nm}, db_{nm}"); b(f"color {col(DIST)}, dist_{nm}"); b(f"set dash_gap, 0.25, dist_{nm}"); b(f"set dash_radius, 0.07, dist_{nm}")
            b(f"hide everything, da_{nm} or db_{nm}"); dist_objs += [f"dist_{nm}", f"da_{nm}", f"db_{nm}"]
        if p.get("polar"):
            b(f"distance polar_{nm}, {nm} and (elem N+O), prot and (elem N+O), 3.5, mode=2"); b(f"color {col(DIST)}, polar_{nm}"); b(f"hide labels, polar_{nm}"); polar_objs.append(f"polar_{nm}")
    if panel.get("reference"):
        files["reference_ligand.pdb"] = panel["reference"]["data"]; b("load reference_ligand.pdb, reference_ligand"); b("show sticks, reference_ligand"); b(f"color {col('#5be37d')}, reference_ligand and elem C"); b("set stick_radius, 0.12, reference_ligand")

    # assemble: colours first
    add(f"set_color c_bg, {_rgb(BG_DARK if dark else BG_LIGHT)}")
    for n, h in colors.items(): add(f"set_color {n}, {_rgb(h)}")
    L.extend(body)

    all_pk = " or ".join(pk_objs) if pk_objs else "prot"
    all_c = " or ".join(cmp_objs)
    fade = lambda t: f"set transparency, {t}, prot"
    add("orient prot")
    # scene 1: overview
    add(fade(0.0)); add("set cartoon_transparency, 0"); add("zoom prot"); add("scene 1_overview, store, view=1, color=1, active=1, rep=1, frame=0")
    # scene 2: pockets highlighted – protein faded (as the explorer's Pocket highlight)
    add(fade(0.6)); add(f"zoom ({all_pk}), 6"); add("scene 2_pockets, store")
    # scene 3: compounds
    if cmp_objs:
        add(fade(0.6)); add(f"zoom ({all_c}), 8"); add("scene 3_compounds, store")
    # scene 4: covalent site
    if tgt:
        add(fade(0.6)); add(f"zoom (target_res{(' or ' + all_c) if cmp_objs else ''}), 6"); add("scene 4_covalent_site, store")
    add(fade(0.0)); add("zoom prot"); add("scene 1_overview, recall")
    add("set_key F1, cmd.scene, 1_overview, recall"); add("set_key F2, cmd.scene, 2_pockets, recall")
    if cmp_objs: add("set_key F3, cmd.scene, 3_compounds, recall")
    if tgt: add("set_key F4, cmd.scene, 4_covalent_site, recall")
    add("# Predictions only: pockets and docked poses are computational and are not evidence of binding, engagement or potency.")
    files["pscreen_scene.pml"] = "\n".join(L) + "\n"
    files["README.txt"] = ("pscreen → PyMOL export\n\n1. Unzip this folder.\n2. In PyMOL: File ▸ Working Directory ▸ Change… (choose this folder), then  @pscreen_scene.pml\n"
                           "   (or: pymol pscreen_scene.pml from a terminal inside this folder).\n3. Scenes: F1 overview · F2 pockets · F3 compounds · F4 covalent site (if a target residue was set).\n\n"
                           "Colours follow the pscreen 3D explorer: pocket colours per prediction tool (P2Rank orange, PocketMiner green, Lacuna blue), target residue yellow,\n"
                           "warhead atoms light orange, reactive atom red, warhead→target distance yellow dashes.\n"
                           "Docked poses and pockets are computational predictions, not evidence of binding, engagement or potency.\n")
    return files


def bundle(panel: dict, dark: bool = True) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in build(panel, dark).items(): z.writestr("pscreen_pymol/" + name, text)
    return buf.getvalue()
