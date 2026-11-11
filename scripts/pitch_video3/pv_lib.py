import time
from playwright.sync_api import sync_playwright
BASE="http://localhost:8501/"
def open_page(p, url, w=1512, h=860):
    b=p.chromium.launch(channel="chrome",headless=True); pg=b.new_page(viewport={"width":w,"height":h}); pg.goto(BASE+url); pg.wait_for_timeout(7000); return b,pg
def xframe(pg):
    for f in pg.frames:
        try:
            if f.evaluate("!!window.__P"): return f
        except Exception: pass
def pill(pg, name):
    pg.get_by_text(name, exact=True).last.click(); pg.wait_for_timeout(3500)
def tick_row(pg, y=460, x=835):
    pg.mouse.click(x,y); pg.wait_for_timeout(3500)
def highlight_target(pg, on=True):
    f=xframe(pg); f.evaluate("""(on)=>{const i=document.querySelector('#hl input[data-k="target"]'); if(i && i.checked!==on){i.checked=on; i.dispatchEvent(new Event('change'));}}""", on); pg.wait_for_timeout(2500)
def focus_camera(pg, xyz, dist=38):
    f=xframe(pg); f.evaluate("""([t,d])=>{const cam=window.__P[0].viewer.plugin.canvas3d.camera; const s=cam.getSnapshot(); const v=[0,1,2].map(i=>s.position[i]-s.target[i]); const l=Math.hypot(...v); cam.setState({target:t, position:[0,1,2].map(i=>t[i]+v[i]/l*d)},0);}""",[xyz,dist]); pg.wait_for_timeout(1500)
