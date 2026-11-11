import sys; sys.path.insert(0,'.')
from groups import *
def ensure_pose(pg,tries=3):
    for _ in range(tries):
        if "to CYS727" in ftext(pg): return
        pg.mouse.click(835,460)
        for _ in range(16):
            pg.wait_for_timeout(250)
            if "to CYS727" in ftext(pg): return
def caps(pg,items,dur):
    """schedule captions (fraction of segment, text) inside the explorer iframe"""
    xframe(pg).evaluate("([it,d])=>{it.forEach(([f,t])=>setTimeout(()=>window.__say('',t),f*d*1000))}",[items,dur])
def g0a(p):
    r=Rec(p,f"?page=Results&run={RUN}","g0a"); pg=r.pg; f=xframe(pg)
    f.get_by_role("button",name="Pockets").first.click(); f.locator("#bNone").click(); pg.keyboard.press("Escape"); pg.wait_for_timeout(1500)
    fs_on(pg); f.evaluate(CAP_JS); f.evaluate("document.querySelector('#bReset').click()"); pg.wait_for_timeout(1500); spin(pg,0.22); r.t0=time.time()-r.c0
    def o1(pg,rem):
        fade(pg,0); d=rem(); caps(pg,[(0.0,"A closed protein structure."),(0.27,"A list of molecules."),(0.45,"One cysteine that matters."),(0.74,"No obvious pocket to dock into.")],d-1.2)
        hold(pg,rem,1.3); fade(pg,1)
    r.seg("o1",o1); r.finish()
def g0b(p):
    cen=pose_centroid(); r=Rec(p,f"?page=Results&run={RUN}","g0b"); pg=r.pg
    ensure_pose(pg); f=xframe(pg); f.get_by_role("button",name="Pockets").first.click(); f.locator("#bNone").click(); pg.keyboard.press("Escape"); pg.wait_for_timeout(800)
    fs_on(pg); xframe(pg).evaluate(CAP_JS); xframe(pg).evaluate("document.querySelector('#bReset').click()"); pg.wait_for_timeout(1500); spin(pg,0.22); r.t0=time.time()-r.c0
    def o2(pg,rem):
        d=rem(); fade(pg,0)
        caps(pg,[(0.0,"Candidate pockets, from independent models."),(0.30,"Covalent poses that reach the cysteine."),(0.55,"Ranked, documented, repeatable."),(0.78,"pscreen: a workbench for structure-based screening.")],d-1.2)
        pg.wait_for_timeout(1200); xframe(pg).evaluate("document.querySelector('#bAll').click()")
        pg.wait_for_timeout(int(d*0.30*1000)-1200)
        spin(pg,0); xframe(pg).evaluate("""()=>{for(const k of ['ligand','target']){const i=document.querySelector(`#hl input[data-k="${k}"]`); if(i&&!i.checked){i.checked=true;i.dispatchEvent(new Event('change'))}}}"""); fly(pg,cen,38,3500)
        pg.wait_for_timeout(int(d*0.25*1000)); fly(pg,cen,85,3500); pg.wait_for_timeout(3600); spin(pg,0.2); hold(pg,rem,1.2); fade(pg,1)
    r.seg("o2",o2); r.finish()
def g1(p):
    r=Rec(p,"","g1")
    def intro(pg,rem):
        t=rem()-0.3; glide(pg,[(183,191),(183,313),(183,558),(183,647)],t*0.5,steps=25)
        per=t*0.5/5
        for x,y in [(390,165),(390,300),(580,190),(770,175),(960,170)]: mv(pg,x,y,22,max(0,int(per*1000)-450))
    def presets(pg,rem):
        t=rem(); mv(pg,743,88,30,400); pg.mouse.click(743,88); pg.wait_for_timeout(900)
        for o in ["Blank (docking only)","P2Rank + PocketMiner + docking + ADMET","Covalent: P2Rank pocket + covalent docking (quick, no ADMET)"]:
            try: bb=pg.get_by_text(o,exact=True).last.bounding_box(); mv(pg,bb["x"]+80,bb["y"]+8,16,500)
            except Exception: pass
        pg.get_by_text("Covalent: cryptic pockets (PocketMiner + Lacuna) + ADMET",exact=True).last.click(); pg.wait_for_timeout(2000); hold(pg,rem)
    def flow(pg,rem):
        t=rem(); glide(pg,[(480,170),(560,210),(680,260),(790,330),(900,200),(1000,170)],t-1.6,steps=30)
        wfr(pg).locator("#bClr").click(); hold(pg,rem,0)
    def build(pg,rem):
        drag_node(pg,"p2rank",(330,420))
        car(pg,"pocketminer"); wf=wfr(pg); mv(pg,205,433,18,350); mv(pg,205,464,12,500)
        wf.locator('.pp[data-t="pocketminer"] .addb').click(); pg.wait_for_timeout(700); car(pg,"pocketminer")
        drag_node(pg,"lacuna",(500,520)); mv(pg,160,617,20,800); mv(pg,160,658,12,800); hold(pg,rem)
    def cov(pg,rem):
        car(pg,"target"); wf=wfr(pg); inp=wf.locator('.pp[data-t="target"] input').first
        bb=inp.bounding_box(); mv(pg,bb["x"]+40,bb["y"]+10,20,250); inp.click(); inp.fill("A:CYS727"); pg.wait_for_timeout(800)
        wf.locator('.pp[data-t="target"] .addb').click(); pg.wait_for_timeout(700); car(pg,"target")
        drag_node(pg,"gnina_cov",(650,560)); hold(pg,rem)
    def admet(pg,rem): drag_node(pg,"admet",(620,520)); wfr(pg).locator("#bFit").click(); hold(pg,rem)
    def inputs(pg,rem):
        mv(pg,1270,140,25,200); pg.get_by_placeholder("e.g. WRN apo · covalent test").fill("WRN · Cys727 · cryptic + covalent"); pg.keyboard.press("Tab"); pg.wait_for_timeout(900)
        mv(pg,1270,210,20,200); pg.locator('input[type=file]').nth(0).set_input_files(ROOT+"data/wrn/9S1A_chainA_apo.pdb"); pg.wait_for_timeout(1800)
        mv(pg,1270,290,20,200); pg.locator('input[type=file]').nth(1).set_input_files(ROOT+"data/wrn/wrn_lib_vvd3.csv"); pg.wait_for_timeout(1800)
        mv(pg,1270,530,20,400); hold(pg,rem,0.8); mv(pg,1274,795,20,400); hold(pg,rem,0)
    for k,f in [("n_intro",intro),("n_presets",presets),("n_flow",flow),("n_build",build),("n_cov",cov),("n_admet",admet),("n_inputs",inputs)]: r.seg(k,f)
    r.finish()
def g2(p):
    r=Rec(p,"?page=Runs","g2")
    def runs(pg,rem): glide(pg,[(300,262),(300,394),(300,437),(300,481),(300,525),(300,569),(1340,150)],rem()-0.3,steps=25)
    r.seg("r_runs",runs); r.finish()
def g3(p):
    cen=pose_centroid(); r=Rec(p,f"?page=Results&run={RUN}","g3")
    def opn(pg,rem): glide(pg,[(430,560),(1100,460),(430,600)],rem()-0.3,steps=25)
    def rotz(pg,rem):
        t=rem(); mv(pg,300,600,15,150); drag_rot(pg,rem,t*0.55); mv(pg,430,600,12,100); n=int((rem()-0.4)/0.4)
        for i in range(n): pg.mouse.wheel(0,-180 if i<n//2 else 180); pg.wait_for_timeout(400)
    def menus(pg,rem):
        t=rem()
        for name in ("Pockets","Residues","View"):
            xframe(pg).get_by_role("button",name=name).first.click(); pg.wait_for_timeout(int((t-0.8)/3*1000)); pg.keyboard.press("Escape"); pg.wait_for_timeout(250)
    def pockets(pg,rem):
        mv(pg,468,240,20,300); pg.get_by_text("Lacuna",exact=True).last.click(); pg.wait_for_timeout(3000); mv(pg,300,600,15,100); drag_rot(pg,rem,rem()-0.3,200)
    def conf(pg,rem):
        sb=pg.locator('[data-testid="stSelectbox"]').filter(has_text="conformer").first
        bx=sb.bounding_box(); mv(pg,bx["x"]+300,bx["y"]+bx["height"]/2,20,200)
        for pick in (3,"docked"):
            sb.click(); pg.wait_for_timeout(700); o=pg.locator('[role="option"]')
            (o.filter(has_text="docked").first if pick=="docked" else o.nth(pick)).click(); pg.wait_for_timeout(1700)
    def pose(pg,rem):
        ensure_pose(pg); focus_camera(pg,cen,34); highlight_target(pg,True); mv(pg,300,600,15,100); drag_rot(pg,rem,rem()-0.3,120)
    def table(pg,rem):
        mv(pg,915,309,20,300); pg.get_by_text("Non-cov.",exact=True).last.click(); pg.wait_for_timeout(max(1500,int(rem()*500)))
        pg.get_by_text("Covalent",exact=True).last.click(); pg.wait_for_timeout(1000); ensure_pose(pg); hold(pg,rem,0)
    def out(pg,rem):
        focus_camera(pg,cen,34); t=rem(); glide(pg,[(860,830),(1006,830),(430,560)],min(6,t*0.4)); spin(pg,0.15)
        pg.evaluate("""()=>{const d=document.createElement('div');d.id='__fade';d.style.cssText='position:fixed;inset:0;background:#05051a;z-index:2147483646;opacity:0;transition:opacity 1.6s;pointer-events:none';document.body.appendChild(d)}""")
        hold(pg,rem,1.8); pg.evaluate("document.getElementById('__fade').style.opacity=1"); pg.wait_for_timeout(1700)
    for k,f in [("r_open",opn),("r_rotz",rotz),("r_menus",menus),("r_pockets",pockets),("r_conf",conf),("r_pose",pose),("r_table",table),("r_out",out)]: r.seg(k,f)
    r.finish()
if __name__=="__main__":
    with sync_playwright() as p:
        for g in sys.argv[1:]: globals()[g](p)
