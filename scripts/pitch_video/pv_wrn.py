import sys; sys.path.insert(0,'.')
from pv_lib import *
import numpy as np, pandas as pd
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit import Chem
res=pd.read_csv("/Users/umairkhan/proteinscreen/runs/WRN_Cys727_cryptic_covalent/results_all_poses.csv")
top=res[(res.docking_mode=="covalent")&(res.covalent_geometry=="covalent_bond_modelled")&(res.compound_id=="VVD-133214")].sort_values("gnina_cnn_vs",ascending=False).iloc[0]
m=[x for x in Chem.SDMolSupplier(top.pose_sdf,removeHs=True,sanitize=False)][int(top.pose_idx)]; cen=m.GetConformer().GetPositions().mean(0).tolist()
with sync_playwright() as p:
    b,pg=open_page(p,"?page=Results&run=WRN_Cys727_cryptic_covalent")
    pill(pg,"Lacuna")                                  # PocketMiner + Lacuna together
    pg.screenshot(path="pv/wrn_overlap.png")
    tick_row(pg); pass
    focus_camera(pg,cen,34)
    pg.screenshot(path="pv/wrn_pose.png")
    b.close()
