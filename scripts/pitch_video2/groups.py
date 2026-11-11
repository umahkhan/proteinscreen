import sys; sys.path.insert(0,'.')
from build import *
RUN="WRN_Cys727_cryptic_covalent"
def g1(p):
    r=Rec(p,"","g1")
    def intro(pg,rem):
        glide(pg,[(357,150),(357,243),(742,150),(1000,160)],rem()-0.5)
    def workflow(pg,rem):
        pg.mouse.move(743,88,steps=30); pg.wait_for_timeout(600); pg.mouse.click(743,88); pg.wait_for_timeout(1800)
        pg.get_by_text("Covalent: cryptic pockets (PocketMiner + Lacuna) + ADMET",exact=True).last.click(); pg.wait_for_timeout(3500)
        glide(pg,[(487,150),(487,243),(613,150),(742,245),(1000,160),(1200,300)],rem()-0.5)
    r.seg("intro",intro); r.seg("workflow",workflow); r.finish()
def g2(p):
    r=Rec(p,"?page=Runs","g2")
    def runs(pg,rem):
        glide(pg,[(300,262),(300,394),(300,437),(300,481),(300,525),(300,569),(1340,150)],rem()-0.5,steps=25)
    r.seg("runs",runs); r.finish()

import pandas as pd
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
def pose_centroid():
    res=pd.read_csv(f"/Users/umairkhan/proteinscreen/runs/{RUN}/results_all_poses.csv")
    top=res[(res.docking_mode=="covalent")&(res.covalent_geometry=="covalent_bond_modelled")&(res.compound_id=="VVD-133214")].sort_values("gnina_cnn_vs",ascending=False).iloc[0]
    m=[x for x in Chem.SDMolSupplier(top.pose_sdf,removeHs=True,sanitize=False)][int(top.pose_idx)]; return m.GetConformer().GetPositions().mean(0).tolist()
def g3(p):
    cen=pose_centroid()
    r=Rec(p,f"?page=Results&run={RUN}","g3")
    def explore(pg,rem):
        spin(pg,0.4); glide(pg,[(430,560),(600,500),(300,640)],rem()-0.8,steps=40); spin(pg,0)
    def pockets(pg,rem):
        pg.mouse.move(468,240,steps=30); pg.wait_for_timeout(500)
        pg.get_by_text("Lacuna",exact=True).last.click(); pg.wait_for_timeout(3500)
        spin(pg,0.3); glide(pg,[(430,560),(250,500)],rem()-0.5,steps=40); spin(pg,0)
    def conformers(pg,rem):
        sb=pg.locator('[data-testid="stSelectbox"]').filter(has_text="conformer").first
        box=sb.bounding_box(); pg.mouse.move(box["x"]+300,box["y"]+box["height"]/2,steps=30); pg.wait_for_timeout(400)
        for pick in (3,6,"docked"):
            sb.click(); pg.wait_for_timeout(900)
            opts=pg.locator('[role="option"]')
            (opts.filter(has_text="docked").first if pick=="docked" else opts.nth(pick)).click(); pg.wait_for_timeout(2400)
    def pose(pg,rem):
        pg.wait_for_timeout(1500)
        for _ in range(3):
            if "to CYS727" in xframe(pg).evaluate("document.body.innerText"): break
            tick_row(pg,y=460,x=835)
        focus_camera(pg,cen,34); highlight_target(pg,True)
        spin(pg,0.15); glide(pg,[(430,520),(520,480)],rem()-0.5,steps=40); spin(pg,0)
    def rank(pg,rem):
        pg.mouse.move(915,309,steps=30); pg.wait_for_timeout(500)
        pg.get_by_text("Non-cov.",exact=True).last.click(); pg.wait_for_timeout(max(2500,int(rem()*500)))
        pg.get_by_text("Covalent",exact=True).last.click(); pg.wait_for_timeout(1500)
        glide(pg,[(1100,460),(1100,560)],rem()-0.3)
    def export(pg,rem):
        pg.wait_for_timeout(500)
        for _ in range(3):
            if "to CYS727" in xframe(pg).evaluate("document.body.innerText"): break
            tick_row(pg,y=460,x=835)
        focus_camera(pg,cen,34); highlight_target(pg,True)
        glide(pg,[(1250,360),(860,830),(1006,830)],max(3,rem()-0.5))
    r.seg("explore",explore); r.seg("pockets",pockets); r.seg("conformers",conformers); r.seg("pose",pose,lead=7); r.seg("rank",rank); r.seg("export",export,lead=7); r.finish()
if __name__=="__main__":
    with sync_playwright() as p:
        for g in sys.argv[1:]: globals()[g](p)
