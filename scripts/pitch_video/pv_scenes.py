import sys; sys.path.insert(0,'.'); sys.path.insert(0,'/Users/umairkhan/proteinscreen')
from pv_lib import *
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
from pathlib import Path
import pandas as pd
from pscreen import api
from pscreen.runio import load_run
def top_cov_centroid(run):
    R=load_run(Path("/Users/umairkhan/proteinscreen/runs/"+run)); rk=api.ranked_compounds(R["results"]); r=rk[rk.docking_mode=="covalent"].iloc[0]
    res=R["results"]; c=res[(res.compound_id==r.compound_id)&(res.docking_mode=="covalent")&(res.pocket_id==r.pocket_id)].sort_values("gnina_cnn_vs",ascending=False).iloc[0]
    m=[x for x in Chem.SDMolSupplier(c.pose_sdf,removeHs=True,sanitize=False)][int(c.pose_idx)]; return m.GetConformer().GetPositions().mean(0).tolist()
SC={"WRN":("WRN_Cys727_cryptic_covalent",["Lacuna"],"wrn"),"PI3K":("PI3Ka_Cys242_cryptic_covalent",["Lacuna"],"pi3k"),"KRAS":("KRAS_G12C_Cys12_all_covalent",["PocketMiner","Lacuna"],"kras")}
which=sys.argv[1:] or list(SC)
with sync_playwright() as p:
    for k in which:
        run,pills,tag=SC[k]; cen=top_cov_centroid(run)
        b,pg=open_page(p,f"?page=Results&run={run}")
        pg.screenshot(path=f"pv/{tag}_single.png")
        for n in pills: pill(pg,n)
        pg.screenshot(path=f"pv/{tag}_overlap.png")
        tick_row(pg); focus_camera(pg,cen,34); pg.screenshot(path=f"pv/{tag}_pose.png")
        b.close()
