import sys; sys.path.insert(0,'.')
from build import *
import pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
ROOT="/Users/umairkhan/proteinscreen/"; RUN="WRN_Cys727_cryptic_covalent"
def ftext(pg):
    try:
        f=xframe(pg); return f.evaluate("document.body.innerText") if f else ""
    except Exception: return ""
def wfr(pg): return [f for f in pg.frames if f.query_selector("#pal")][0]
def mv(pg,x,y,steps=30,wait=400): pg.mouse.move(x,y,steps=steps); pg.wait_for_timeout(wait)
def hold(pg,rem,extra=0.3):
    t=rem()-extra
    if t>0: pg.wait_for_timeout(int(t*1000))
def drag_node(pg,t,pos):
    wf=wfr(pg); src=wf.locator(f'.pi[data-t="{t}"]'); src.scroll_into_view_if_needed(); pg.wait_for_timeout(300)
    bb=src.bounding_box(); mv(pg,bb["x"]+60,bb["y"]+bb["height"]/2,25,300)
    src.drag_to(wf.locator("#cv"),target_position={"x":pos[0],"y":pos[1]}); pg.wait_for_timeout(900)
def car(pg,t): wfr(pg).locator(f'.pi[data-t="{t}"] .car').click(); pg.wait_for_timeout(700)

def g1(p):
    r=Rec(p,"","g1")
    def intro(pg,rem): glide(pg,[(467,33),(571,33),(669,33),(300,120)],rem()-0.3,steps=30)
    def palette(pg,rem): glide(pg,[(183,191),(183,313),(183,558),(183,647),(600,400)],rem()-0.3,steps=30)
    def essentials(pg,rem):
        t=rem(); per=(t-0.5)/5
        for x,y in [(390,165),(390,300),(580,190),(770,175),(960,170)]:
            mv(pg,x,y,25,max(0,int(per*1000)-500))
    def presets(pg,rem):
        t=rem(); mv(pg,743,88,30,500); pg.mouse.click(743,88); pg.wait_for_timeout(1200)
        opts=["Blank (docking only)","P2Rank + docking (no ADMET)","P2Rank + PocketMiner + docking + ADMET","Covalent: P2Rank pocket + covalent docking (quick, no ADMET)"]
        for o in opts:
            try: bb=pg.get_by_text(o,exact=True).last.bounding_box(); mv(pg,bb["x"]+80,bb["y"]+8,18,max(200,int((t*0.45-2)/4*1000)))
            except Exception: pass
        pg.get_by_text("Covalent: cryptic pockets (PocketMiner + Lacuna) + ADMET",exact=True).last.click(); pg.wait_for_timeout(2500)
        hold(pg,rem)
    def flow(pg,rem):
        t=rem(); glide(pg,[(480,170),(560,210),(680,260),(790,330),(900,200),(1000,170)],t-1.5,steps=35)
        wfr(pg).locator("#bClr").click(); hold(pg,rem,0)
    def site(pg,rem):
        mv(pg,183,313,25,600); mv(pg,229,350,20,700)    # tooltip on P2Rank '?'
        hold(pg,rem,3.0); drag_node(pg,"p2rank",(330,420)); hold(pg,rem)
    def pocketminer(pg,rem):
        car(pg,"pocketminer"); wf=wfr(pg)
        mv(pg,205,433,25,600); mv(pg,205,464,15,900)   # params
        wf.locator('.pp[data-t="pocketminer"] .addb').click(); pg.wait_for_timeout(1000); car(pg,"pocketminer"); hold(pg,rem)
    def lacuna(pg,rem): drag_node(pg,"lacuna",(500,520)); hold(pg,rem)
    def known(pg,rem):
        mv(pg,160,617,25,1500); mv(pg,160,658,15,1500); hold(pg,rem)
    def target(pg,rem):
        car(pg,"target"); wf=wfr(pg); inp=wf.locator('.pp[data-t="target"] input').first
        bb=inp.bounding_box(); mv(pg,bb["x"]+40,bb["y"]+10,25,300); inp.click(); inp.fill("A:CYS727"); pg.wait_for_timeout(1200)
        wf.locator('.pp[data-t="target"] .addb').click(); pg.wait_for_timeout(1000); car(pg,"target"); hold(pg,rem)
    def covalent(pg,rem): drag_node(pg,"gnina_cov",(650,560)); hold(pg,rem)
    def admet(pg,rem): drag_node(pg,"admet",(620,520)); hold(pg,rem)
    def review(pg,rem):
        wfr(pg).locator("#bFit").click(); pg.wait_for_timeout(800)
        glide(pg,[(400,200),(700,250),(900,450),(500,500),(1000,700)],rem()-0.3,steps=30)
    def inputs(pg,rem):
        mv(pg,1270,140,30,300); pg.get_by_placeholder("e.g. WRN apo · covalent test").fill("WRN · Cys727 · cryptic + covalent"); pg.keyboard.press("Tab"); pg.wait_for_timeout(1200)
        mv(pg,1270,210,25,300); pg.locator('input[type=file]').nth(0).set_input_files(ROOT+"data/wrn/9S1A_chainA_apo.pdb"); pg.wait_for_timeout(2500)
        mv(pg,1270,290,25,300); pg.locator('input[type=file]').nth(1).set_input_files(ROOT+"data/wrn/wrn_lib_vvd3.csv"); pg.wait_for_timeout(2500)
        mv(pg,1270,530,25,500); hold(pg,rem,1.2); mv(pg,1274,795,25,500); hold(pg,rem,0)
    for k,f in [("intro",intro),("palette",palette),("essentials",essentials),("presets",presets),("flow",flow),("site",site),("pocketminer",pocketminer),("lacuna",lacuna),("known",known),("target",target),("covalent",covalent),("admet",admet),("review",review),("inputs",inputs)]:
        r.seg(k,f)
    r.finish()
def g2(p):
    r=Rec(p,"?page=Runs","g2")
    def runs(pg,rem): glide(pg,[(300,262),(300,394),(300,437),(300,481),(300,525),(300,569),(1340,150)],rem()-0.3,steps=25)
    r.seg("runs",runs); r.finish()
def pose_centroid():
    res=pd.read_csv(f"{ROOT}runs/{RUN}/results_all_poses.csv")
    top=res[(res.docking_mode=="covalent")&(res.covalent_geometry=="covalent_bond_modelled")&(res.compound_id=="VVD-133214")].sort_values("gnina_cnn_vs",ascending=False).iloc[0]
    m=[x for x in Chem.SDMolSupplier(top.pose_sdf,removeHs=True,sanitize=False)][int(top.pose_idx)]; return m.GetConformer().GetPositions().mean(0).tolist()
def drag_rot(pg,rem,total,dx=260):
    t_end=time.time()+total; sgn=1; x0=430; y=600
    while time.time()<t_end-0.2:
        pg.mouse.move(x0-sgn*dx/2,y); pg.mouse.down(); pg.mouse.move(x0+sgn*dx/2,y+20,steps=40); pg.mouse.up(); sgn=-sgn; pg.wait_for_timeout(150)

def fs_on(pg):
    pg.evaluate("""()=>{if(!document.getElementById('__fade')){const d=document.createElement('div');d.id='__fade';d.style.cssText='position:fixed;inset:0;background:#05051a;z-index:2147483646;opacity:1;transition:opacity .9s;pointer-events:none';document.body.appendChild(d);}}""")
    pg.evaluate("""()=>{if(!document.getElementById('__bg')){const d=document.createElement('div');d.id='__bg';d.style.cssText='position:fixed;inset:0;background:#05051a;z-index:2147483000';document.body.appendChild(d);}}""")
    xframe(pg).evaluate("()=>{if(!document.getElementById('__hb')){const s=document.createElement('style');s.id='__hb';s.textContent='#bar{display:none!important}';document.head.appendChild(s);}}")
    xframe(pg).frame_element().evaluate("e=>{e.dataset.o=e.style.cssText;e.style.cssText+=';position:fixed!important;inset:0!important;width:100vw!important;height:100vh!important;z-index:2147483640!important;background:#05051a'}")
def fade(pg,to): pg.evaluate("(o)=>{document.getElementById('__fade').style.opacity=o}",to); pg.wait_for_timeout(1000)
CAP_JS="""(()=>{if(document.getElementById('__cap'))return;const st=document.createElement('style');st.textContent=`
#__cap{position:fixed;inset:0;pointer-events:none;z-index:99999;background:radial-gradient(ellipse at 50% 45%,rgba(5,5,26,0) 35%,rgba(5,5,26,.78) 100%)}
#__cap .box{position:absolute;left:7vw;bottom:11vh;max-width:70vw}
#__cap .k{font:800 22px 'Helvetica Neue',Arial;letter-spacing:.28em;color:#3fd0c0;opacity:0;transform:translateY(14px);transition:all .8s}
#__cap .l{font:800 62px/1.08 'Helvetica Neue',Arial;color:#fff;text-shadow:0 4px 30px rgba(0,0,0,.8);opacity:0;transform:translateY(22px);filter:blur(8px);transition:all .9s cubic-bezier(.2,.7,.2,1);margin-top:8px}
#__cap .bar{height:4px;width:0;background:#3fd0c0;margin-top:18px;transition:width 1.2s}
#__cap .on{opacity:1!important;transform:none!important;filter:none!important}`;document.head.appendChild(st);
const c=document.createElement('div');c.id='__cap';c.innerHTML='<div class="box"><div class="k"></div><div class="l"></div><div class="bar"></div></div>';document.body.appendChild(c);
window.__say=(k,l)=>{const K=c.querySelector('.k'),L=c.querySelector('.l'),B=c.querySelector('.bar');[K,L].forEach(e=>e.classList.remove('on'));B.style.width=0;
 setTimeout(()=>{K.textContent=k;L.textContent=l;K.classList.add('on');L.classList.add('on');B.style.width='140px';},650);};})()"""
def say(pg,k,l): xframe(pg).evaluate("([k,l])=>window.__say(k,l)",[k,l])
def fly(pg,t,d,ms=3500):
    """Fly the camera to look at t from outside the protein (direction = protein centre -> t), ending d A away."""
    xframe(pg).evaluate("""([t,d,ms])=>{const cam=window.__P[0].viewer.plugin.canvas3d.camera;const s=cam.getSnapshot();
 if(!window.__C)window.__C=s.target.slice();const C=window.__C;
 const nrm=v=>{const l=Math.hypot(...v);return v.map(x=>x/l)};
 const u1=nrm([0,1,2].map(i=>t[i]-C[i]));const v0=[0,1,2].map(i=>s.position[i]-s.target[i]);const D0=Math.hypot(...v0),u0=nrm(v0),T0=s.target.slice(),t0=performance.now();
 const f=n=>{let a=Math.min(1,(n-t0)/ms);a=a*a*(3-2*a);const tt=[0,1,2].map(i=>T0[i]+(t[i]-T0[i])*a),u=nrm([0,1,2].map(i=>u0[i]+(u1[i]-u0[i])*a)),dd=D0+(d-D0)*a;
  cam.setState({target:tt,position:[0,1,2].map(i=>tt[i]+u[i]*dd)},0);if(a<1)requestAnimationFrame(f)};requestAnimationFrame(f)}""",[t,d,ms])
def g0(p):
    cen=pose_centroid(); r=Rec(p,f"?page=Results&run={RUN}","g0")
    pg=r.pg; f=xframe(pg)
    # stage the BEFORE look: no pockets, no pose
    f.locator("#bNone").evaluate("e=>e.click()") if False else None
    f.get_by_role("button",name="Pockets").first.click(); f.locator("#bNone").click(); pg.keyboard.press("Escape"); pg.wait_for_timeout(1500)
    fs_on(pg); f.evaluate(CAP_JS); f.evaluate("document.querySelector('#bReset').click()"); pg.wait_for_timeout(1500); spin(pg,0.22); r.t0=time.time()-r.c0
    def b1(pg,rem): fade(pg,0); say(pg,"BEFORE","A closed protein structure."); hold(pg,rem,0)
    def b2(pg,rem): say(pg,"BEFORE","A list of molecules."); hold(pg,rem,0)
    def b3(pg,rem): say(pg,"BEFORE","And one cysteine you believe matters."); hold(pg,rem,0)
    def b4(pg,rem): say(pg,"BEFORE","But no obvious pocket to dock into."); hold(pg,rem,1.2); fade(pg,1); spin(pg,0)
    def setup_after(pg):
        pg.evaluate("""()=>{document.querySelectorAll('iframe').forEach(e=>{if(e.dataset.o!==undefined){e.style.cssText=e.dataset.o;delete e.dataset.o}});const b=document.getElementById('__bg');if(b)b.remove()}""")
        pg.wait_for_timeout(1200)
        pg.mouse.click(835,460)
        for _ in range(25):
            pg.wait_for_timeout(300)
            if "to CYS727" in ftext(pg): break
        pg.wait_for_timeout(600); f2=xframe(pg); f2.get_by_role("button",name="Pockets").first.click(); f2.locator("#bNone").click(); pg.keyboard.press("Escape"); pg.wait_for_timeout(500)
        fs_on(pg); xframe(pg).evaluate("document.querySelector('#bReset').click()"); pg.wait_for_timeout(1500); spin(pg,0.22)
    def a1(pg,rem):
        fade(pg,0); say(pg,"AFTER","Candidate pockets, found by independent models.")
        pg.wait_for_timeout(1500); xframe(pg).evaluate("document.querySelector('#bAll').click()") ; hold(pg,rem,0)
    def a2(pg,rem):
        say(pg,"AFTER","A covalent pose, reaching the cysteine."); spin(pg,0)
        xframe(pg).evaluate("""()=>{for(const k of ['ligand','target']){const i=document.querySelector(`#hl input[data-k="${k}"]`); if(i&&!i.checked){i.checked=true;i.dispatchEvent(new Event('change'))}}}"""); fly(pg,cen,38,3500); hold(pg,rem,0)
    def a3(pg,rem):
        say(pg,"AFTER","Ranked, documented, and reproducible."); spin(pg,0); fly(pg,cen,85,3500); pg.wait_for_timeout(3700); spin(pg,0.2); hold(pg,rem,0)
    def a4(pg,rem):
        say(pg,"pscreen","A workbench for structure-based screening."); hold(pg,rem,1.2); fade(pg,1)
    r.seg("b1",b1); r.seg("b2",b2); r.seg("b3",b3); r.seg("b4",b4)
    setup_after(pg); xframe(pg).evaluate(CAP_JS)
    xframe(pg).evaluate("document.querySelector('#bNone').click()") if False else None
    r.seg("a1",a1); r.seg("a2",a2); r.seg("a3",a3); r.seg("a4",a4); r.finish()
def g3(p):
    cen=pose_centroid(); r=Rec(p,f"?page=Results&run={RUN}","g3")
    def results(pg,rem): glide(pg,[(430,560),(1100,460),(430,600)],rem()-0.3,steps=30)
    def rotate(pg,rem): mv(pg,300,600,20,200); drag_rot(pg,rem,rem()-0.3)
    def zoom(pg,rem):
        mv(pg,430,600,15,200); n=int((rem()-0.5)/0.45)
        for i in range(n): pg.mouse.wheel(0,-180 if i<n//2 else 180); pg.wait_for_timeout(450)
    def menus(pg,rem):
        t=rem()
        for name in ("Pockets","Residues","View"):
            xframe(pg).get_by_role("button",name=name).first.click(); pg.wait_for_timeout(int((t-0.8)/3*1000)); pg.keyboard.press("Escape"); pg.wait_for_timeout(250)
    def pockets(pg,rem):
        mv(pg,468,240,25,400); pg.get_by_text("Lacuna",exact=True).last.click(); pg.wait_for_timeout(3500); mv(pg,300,600,20,100); drag_rot(pg,rem,rem()-0.3,200)
    def conformers(pg,rem):
        sb=pg.locator('[data-testid="stSelectbox"]').filter(has_text="conformer").first
        bx=sb.bounding_box(); mv(pg,bx["x"]+300,bx["y"]+bx["height"]/2,25,300)
        for pick in (3,6,"docked"):
            sb.click(); pg.wait_for_timeout(800); o=pg.locator('[role="option"]')
            (o.filter(has_text="docked").first if pick=="docked" else o.nth(pick)).click(); pg.wait_for_timeout(1900)
    def pose(pg,rem):
        for _ in range(3):
            if "to CYS727" in ftext(pg): break
            tick_row(pg,y=460,x=835)
        focus_camera(pg,cen,34); highlight_target(pg,True); mv(pg,300,600,20,100); drag_rot(pg,rem,rem()-0.3,120)
    def table(pg,rem):
        mv(pg,915,309,25,400); pg.get_by_text("Non-cov.",exact=True).last.click(); pg.wait_for_timeout(2500)
        pg.get_by_text("Covalent",exact=True).last.click(); pg.wait_for_timeout(1500); glide(pg,[(1100,460),(1100,560),(1200,354)],max(1,rem()-0.3))
        for _ in range(3):
            if "to CYS727" in ftext(pg): break
            tick_row(pg,y=460,x=835)
    def outro(pg,rem):
        focus_camera(pg,cen,34); glide(pg,[(860,830),(1006,830)],6); hold(pg,rem,0)
    for k,f in [("results",results),("rotate",rotate),("zoom",zoom),("menus",menus),("pockets",pockets),("conformers",conformers),("pose",pose),("table",table),("outro",outro)]:
        r.seg(k,f)
    r.finish()
if __name__=="__main__":
    with sync_playwright() as p:
        for g in sys.argv[1:]: globals()[g](p)
