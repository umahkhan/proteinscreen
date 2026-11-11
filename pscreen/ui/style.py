"""Visual system: one blue accent, orange for the covalent lane, neutral surfaces (validated dataviz palette)."""
import streamlit as st

BLUE, ORANGE, AQUA, VIOLET, RED, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e34948", "#008300", "#8a8a85"
LANE_COLOR = {"noncovalent": BLUE, "covalent": ORANGE}
LOGO = ('<svg width="34" height="34" viewBox="0 0 34 34"><rect width="34" height="34" rx="9" fill="#14161a"/><path d="M17 6.5l8.2 4.75v9.5L17 25.5l-8.2-4.75v-9.5z" fill="none" stroke="#fff" stroke-width="1.8"/>'
        '<circle cx="17" cy="16" r="3" fill="#2a78d6"/><circle cx="23.2" cy="12.4" r="1.6" fill="#eb6834"/><circle cx="10.8" cy="19.6" r="1.6" fill="#1baf7a"/></svg>')
POCKET_COLORS = [BLUE, ORANGE, AQUA]       # first three categorical slots validate all-pairs; more pockets fold to grey in charts

CSS = f"""
<style>
[data-testid="stToolbar"], #MainMenu, footer, [data-testid="stDecoration"] {{ display:none !important; }}
.ps-hd {{ display:flex; align-items:baseline; gap:14px; margin:.2rem 0 .8rem; }}
.ps-t {{ font-size:1.6rem; font-weight:800; letter-spacing:-.01em; color:#1b1b1a; }}
.ps-s {{ font-size:.88rem; font-weight:600; color:#55544f; }}
.ps-div {{ height:1px; background:#ecebe7; margin:.7rem 0 .55rem; }}
[data-testid="stButtonGroup"] {{ flex-wrap:nowrap !important; gap:0 !important; }}
[data-testid="stButtonGroup"] button {{ white-space:nowrap; padding:0.25rem 0.6rem; }}
.st-key-ctl label p {{ font-size:.7rem !important; font-weight:800 !important; letter-spacing:.07em; text-transform:uppercase; color:#1b1b1a !important; }}
.st-key-ctl label {{ min-height:0 !important; margin-bottom:2px !important; }}
.st-key-ctl [data-baseweb="select"] > div, .st-key-ctl [data-testid="stNumberInputContainer"] {{ background:#e9ecf1 !important; border:1.5px solid #8b93a1 !important; font-weight:650; color:#111 !important; min-height:2.1rem; }}
.st-key-ctl [data-testid="stNumberInput"] input {{ font-weight:700; color:#111; }}
.st-key-ctl [data-testid="stVerticalBlock"] {{ gap:.35rem !important; }}
.st-key-ctl [data-testid="stButtonGroup"] button {{ font-weight:700; border:1.5px solid #8b93a1; }}
.stPills button, [data-testid="stPills"] button, [data-baseweb="select"] div, .stMultiSelect, [data-testid="stPopover"] button {{ font-weight:650 !important; }}
[data-testid="stPills"] button, [data-testid="stBaseButton-pills"], [data-testid="stBaseButton-pillsActive"] {{ border:1.5px solid #8b93a1 !important; font-weight:700 !important; color:#111 !important; }}
.st-key-wfin [data-testid="stFileUploaderDropzone"] {{ padding:.35rem .6rem; flex-direction:row; gap:.5rem; min-height:0; border:1.5px dashed #8b93a1; background:#f3f4f6; }}
.st-key-wfin [data-testid="stFileUploaderDropzoneInstructions"] {{ display:none; }}
.st-key-wfin [data-testid="stFileUploaderDropzone"] button {{ padding:.15rem .6rem; font-weight:700; border:1.5px solid #8b93a1; width:100%; }}
.st-key-wfin [data-testid="stFileUploaderFile"] {{ padding:0 !important; min-height:0 !important; }}
.st-key-wfin [data-testid="stFileUploaderFile"] small {{ display:none; }}
.st-key-wfin [data-testid="stFileUploaderFileName"] {{ font-weight:700; font-size:.78rem; }}
.st-key-wf_run {{ margin-bottom:30px; position:relative; top:-14px; }}
.st-key-wfin label p {{ font-size:.68rem !important; font-weight:800 !important; letter-spacing:.06em; text-transform:uppercase; color:#1b1b1a !important; }}
.st-key-wfin label {{ margin-bottom:1px !important; min-height:0 !important; }}
.st-key-wfin [data-testid="stVerticalBlock"] {{ gap:.45rem !important; }}
.st-key-wfin input {{ font-weight:700; }}
.ps-h {{ font-size:.82rem; font-weight:800; letter-spacing:.07em; text-transform:uppercase; color:#1b1b1a; border-left:3px solid #2a78d6; padding-left:9px; margin:.1rem 0 .45rem; line-height:1.1; }}
.ps-gl {{ font-size:.68rem; font-weight:700; letter-spacing:.08em; text-transform:uppercase; color:#6b6a66; margin:0 0 .25rem; }}
html, body, [class*="css"] {{ font-family: -apple-system, "Segoe UI", Inter, Helvetica, Arial, sans-serif; }}
div[data-testid="stVerticalBlockBorderWrapper"] {{ border-radius:12px; border-color:#e2e1dc; }}
:root {{ --blue:{BLUE}; --orange:{ORANGE}; --ink:#1b1b1a; --muted:#5f5e5a; --line:#e6e5e1; --soft:#f6f6f4; }}
section[data-testid="stSidebar"], [data-testid="stSidebarCollapsedControl"], [data-testid="collapsedControl"] {{ display:none !important; }}
.block-container {{ padding-top: 0.5rem; padding-bottom: 0.2rem; max-width: 1800px; }}
header[data-testid="stHeader"] {{ display:none; }}
h1, h2, h3 {{ letter-spacing: -0.01em; }}
.ps-title {{ font-size: 1.55rem; font-weight: 650; margin: 0; color: var(--ink); }}
.ps-sub {{ color: var(--muted); font-size: .86rem; margin: .1rem 0 .9rem; }}
.ps-cards {{ display:grid; grid-template-columns: repeat(auto-fit,minmax(150px,1fr)); gap:10px; margin: .4rem 0 1rem; }}
.ps-card {{ border:1px solid var(--line); border-radius:10px; padding:10px 14px; background:#fff; }}
.ps-card .v {{ font-size:1.55rem; font-weight:650; color:var(--ink); line-height:1.15; }}
.ps-card .k {{ font-size:.74rem; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); }}
.ps-card .n {{ font-size:.76rem; color:var(--muted); margin-top:2px; }}
.ps-badge {{ display:inline-block; padding:1px 8px; margin:0 4px 4px 0; border-radius:999px; font-size:.74rem; border:1px solid var(--line); background:var(--soft); color:var(--ink); }}
.ps-badge.blue {{ background:#e8f0fb; border-color:#c3d7f3; color:#1b5aa8; }}
.ps-badge.orange {{ background:#fdeee7; border-color:#f6cfbf; color:#a8421c; }}
.ps-badge.green {{ background:#e6f5ee; border-color:#bfe3d2; color:#0f6c4a; }}
.ps-badge.red {{ background:#fdeaea; border-color:#f4c5c5; color:#a12c2b; }}
.ps-badge.grey {{ color:var(--muted); }}
.ps-note {{ border-left:3px solid var(--orange); background:#fff8f4; padding:8px 12px; border-radius:4px; font-size:.84rem; color:#5a3a2c; margin:.4rem 0; }}
.ps-info {{ border-left:3px solid var(--blue); background:#f4f8fd; padding:8px 12px; border-radius:4px; font-size:.84rem; color:#23405f; margin:.4rem 0; }}
.ps-legend td, .ps-legend th {{ font-size:.82rem; padding:3px 10px 3px 0; vertical-align:top; }}
.ps-legend th {{ text-align:left; color:var(--muted); font-weight:600; }}
div[data-testid="stTabs"] button[role="tab"] {{ font-size:.95rem; padding: 6px 14px; }}
[data-testid="stSidebar"] {{ background:#fafaf8; border-right:1px solid var(--line); }}

/* ---- production shell ---- */
[data-testid="stElementContainer"]:has(> iframe[height="0"]), [data-testid="stElementContainer"]:has(iframe[height="0"]):not(:has(iframe[height]:not([height="0"]))) {{ display:none !important; }}
[data-testid="stMain"] [data-testid="stVerticalBlock"] {{ gap:.55rem; }}
.ps-brand {{ display:flex; align-items:center; gap:11px; }}
.ps-brand svg {{ flex:none; }}
.ps-brand-t {{ font-size:1.35rem; font-weight:900; letter-spacing:-.02em; color:#14161a; line-height:1.05; }}
.ps-brand-s {{ font-size:.68rem; font-weight:700; letter-spacing:.07em; text-transform:uppercase; color:#5b6270; }}
.ps-live {{ text-align:right; font-size:1.05rem; font-weight:800; color:#3b414c; }}
.ps-live .dot {{ display:inline-block; width:12px; height:12px; border-radius:50%; background:#2a78d6; margin-right:9px; animation:psp 1.2s ease-in-out infinite; }}
@keyframes psp {{ 50% {{ opacity:.3; transform:scale(1.35); }} }}
.st-key-nav [data-testid="stButtonGroup"] button {{ font-size:.92rem; font-weight:800; padding:.5rem 1.35rem; border:1.5px solid #8b93a1; color:#14161a; }}
.st-key-nav [data-testid="stButtonGroup"] button[aria-checked="true"], .st-key-nav [data-testid="stBaseButton-segmented_controlActive"] {{ background:#14161a !important; color:#fff !important; border-color:#14161a !important; }}
.ps-kick {{ font-size:.7rem; font-weight:900; letter-spacing:.11em; text-transform:uppercase; color:#5b6270; }}
.ps-big {{ font-size:1.55rem; font-weight:900; letter-spacing:-.015em; color:#14161a; line-height:1.15; display:flex; align-items:center; gap:12px; flex-wrap:wrap; }}
.ps-meta {{ font-size:.82rem; font-weight:650; color:#4b5260; margin-top:3px; }}
.ps-meta code, .ps-ctx code {{ background:#eef0f4; border-radius:4px; padding:0 5px; font-size:.78rem; font-weight:700; color:#14161a; }}
.ps-card-h {{ display:flex; justify-content:space-between; align-items:flex-start; gap:20px; border:1.5px solid #8b93a1; border-radius:12px; padding:16px 20px; background:#fff; box-shadow:0 1px 3px rgba(20,22,26,.07); margin-bottom:12px; }}
.ps-ctx {{ border:1.5px solid #14161a; border-left-width:7px; border-radius:10px; padding:7px 16px; background:#f7f9fc; line-height:1.15; margin-bottom:.3rem; }}
.ps-ctx-r {{ display:flex; justify-content:space-between; gap:12px; align-items:baseline; }}
.ps-nodes {{ display:flex; flex-wrap:wrap; align-items:center; gap:4px; margin-top:4px; }}
.ps-node {{ display:inline-block; border:1.5px solid #8b93a1; border-left-width:5px; border-radius:6px; padding:1px 8px; font-size:.72rem; font-weight:800; background:#fff; color:#14161a; }}
.ps-arrow {{ color:#8b93a1; font-weight:900; margin:0 1px; }}
.ps-note {{ font-size:.82rem; color:#5b6270; font-weight:600; padding-top:1.1rem; }}
.ps-ctx {{ margin-bottom:1rem; }}
.ps-ctx2 {{ display:flex; align-items:stretch; justify-content:space-between; gap:24px; padding:9px 18px 10px 16px; }}
.ps-ctx2-l {{ display:flex; flex-direction:column; gap:3px; min-width:0; }}
.ps-ctx2-r {{ display:flex; gap:26px; align-items:center; border-left:1.5px solid #d5d9e0; padding-left:24px; flex:none; }}
.ps-fact {{ display:flex; flex-direction:column; gap:3px; }}
.ps-fact span {{ font-size:.66rem; font-weight:900; letter-spacing:.11em; text-transform:uppercase; color:#5b6270; }}
.ps-fact b, .ps-fact code {{ font-size:.86rem; font-weight:800; color:#14161a; }}
.ps-ctx .ps-big {{ font-size:1.2rem; }} .ps-ctx .ps-meta {{ margin-top:1px; font-size:.76rem; }}
.ps-right {{ text-align:right; }}
.ps-timer {{ font-size:2rem; font-weight:900; font-variant-numeric:tabular-nums; letter-spacing:-.02em; color:#14161a; line-height:1.1; }}
.ps-chip {{ display:inline-block; border-radius:999px; padding:3px 12px; font-size:.74rem; font-weight:900; letter-spacing:.03em; }}
.ps-chip.s {{ padding:2px 10px; font-size:.7rem; }}
table.ps-steps {{ width:100%; border-collapse:separate; border-spacing:0; border:1.5px solid #8b93a1; border-radius:12px; overflow:hidden; background:#fff; font-size:.86rem; }}
table.ps-steps th {{ text-align:left; background:#eef0f4; color:#3b414c; font-size:.68rem; font-weight:900; letter-spacing:.09em; text-transform:uppercase; padding:9px 12px; border-bottom:1.5px solid #8b93a1; }}
table.ps-steps td {{ padding:9px 12px; border-bottom:1px solid #e3e6ec; vertical-align:middle; font-weight:650; color:#14161a; }}
table.ps-steps tr:last-child td {{ border-bottom:none; }}
table.ps-steps td.n {{ width:34px; color:#8b93a1; font-weight:900; text-align:center; }}
table.ps-steps td.tm {{ text-align:right; font-variant-numeric:tabular-nums; font-weight:900; white-space:nowrap; }}
table.ps-steps td.msg {{ color:#4b5260; font-size:.78rem; }}
.ps-bar {{ display:inline-block; vertical-align:middle; width:150px; height:8px; border-radius:5px; background:#e3e6ec; overflow:hidden; }}
.ps-bar div {{ height:8px; border-radius:5px; }}
.ps-cnt {{ margin-left:8px; font-size:.76rem; font-weight:800; color:#3b414c; }}
.ps-th {{ font-size:.66rem; font-weight:900; letter-spacing:.09em; text-transform:uppercase; color:#5b6270; border-bottom:1.5px solid #8b93a1; padding-bottom:5px; }}
.ps-row .ps-meta {{ white-space:nowrap; overflow:hidden; text-overflow:ellipsis; margin:2px 0 0; max-width:100%; }}
.ps-row {{ border-left:4px solid transparent; padding-left:10px; line-height:1.25; }}
.ps-row.sel {{ border-left-color:#2a78d6; background:#f1f6fe; border-radius:0 6px 6px 0; }}
.ps-row b {{ font-size:.92rem; font-weight:850; }}
.ps-cell {{ font-size:.82rem; font-weight:700; color:#25292f; }} .ps-cell.tab {{ font-variant-numeric:tabular-nums; font-weight:800; }}
.ps-empty {{ border:1.5px dashed #8b93a1; border-radius:12px; padding:34px; text-align:center; color:#3b414c; font-size:.95rem; font-weight:650; background:#fafbfc; }}
.ps-err {{ font-size:.8rem; font-weight:800; color:#a12c2b; margin:2px 0; }}
.ps-warn {{ font-size:.8rem; font-weight:800; color:#8a5a00; margin:2px 0; }}
.st-key-wfin {{ border:1.5px solid #8b93a1 !important; border-radius:12px; background:#fff; box-shadow:0 1px 3px rgba(20,22,26,.07); padding:10px 14px 10px !important; box-sizing:border-box; }}
.st-key-wfin [data-testid="stVerticalBlock"] {{ gap:.35rem !important; }}
.st-key-wfin label p {{ font-size:.64rem !important; white-space:normal; line-height:1.2; }}
.st-key-wfin label {{ width:100%; justify-content:space-between; gap:6px; }}
.st-key-wfin [data-testid="stWidgetLabel"] {{ min-height:0; padding-right:2px; }}
.st-key-wfin [data-testid="stFileUploaderDropzone"] button {{ padding:.05rem .5rem; font-size:.8rem; }}
.st-key-wfin [data-testid="stCaptionContainer"], .st-key-wfin .ps-err, .st-key-wfin .ps-warn {{ font-size:.72rem !important; }}
</style>
"""


def inject(): st.markdown(CSS, unsafe_allow_html=True)


def cards(items: list[tuple[str, str, str]]):
    """items: (label, value, note)."""
    h = "".join(f'<div class="ps-card"><div class="k">{k}</div><div class="v">{v}</div><div class="n">{n}</div></div>' for k, v, n in items)
    st.markdown(f'<div class="ps-cards">{h}</div>', unsafe_allow_html=True)


def badge(text: str, kind: str = "grey") -> str:
    return f'<span class="ps-badge {kind}">{text}</span>'


def note(text: str, kind: str = "info"): st.markdown(f'<div class="ps-{kind}">{text}</div>', unsafe_allow_html=True)


SCORE_LEGEND = [
    ("Vina-like score (kcal/mol)", "GNINA minimizedAffinity. More negative = better predicted fit. Tends to favour larger ligands.", "gnina_minimized_affinity_kcal"),
    ("CNN pose score (0–1)", "Probability that the pose resembles a good crystal-like pose. It rates the POSE, not the compound.", "gnina_cnn_pose_score"),
    ("CNN affinity (pK)", "Predicted −log10(affinity); higher = tighter. A model estimate, not a measurement.", "gnina_cnn_affinity_pK"),
    ("CNN VS (pose × affinity)", "CNN pose score × CNN affinity: GNINA's single-number ranking option. Used as the default sort.", "gnina_cnn_vs"),
]


