import sys; sys.path.insert(0,'.')
from pv_lib import *
with sync_playwright() as p:
    b,pg=open_page(p,"?page=Workflow")
    pg.get_by_text("Start from a preset").first.click(); pg.wait_for_timeout(1000)
    pg.screenshot(path="pv/wf_menu.png")
    opts=pg.locator('[role="option"]').all_inner_texts(); print(opts)
    for o in opts:
        if "cryptic" in o.lower() and "covalent" in o.lower(): pg.get_by_role("option", name=o).click(); break
    pg.wait_for_timeout(5000); pg.screenshot(path="pv/wf_preset.png")
    b.close()
    b,pg=open_page(p,"?page=Runs&run=WRN_Cys727_cryptic_covalent"); pg.screenshot(path="pv/runs_wrn.png"); b.close()
