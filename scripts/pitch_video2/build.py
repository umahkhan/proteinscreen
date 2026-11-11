import sys, time, subprocess, shutil, os
sys.path.insert(0,'.')
from pv_lib import *
from pron import speak, np, sf
from script import SEGS
import imageio_ffmpeg
FF=imageio_ffmpeg.get_ffmpeg_exe(); OUT="pv/"; os.makedirs(OUT,exist_ok=True)
CURSOR="""(()=>{const mk=()=>{if(document.getElementById('__cur'))return null;const d=document.createElement('div');d.id='__cur';
d.style.cssText='position:fixed;z-index:2147483647;width:20px;height:20px;border-radius:50%;background:rgba(255,196,0,.55);border:2px solid #fff;box-shadow:0 0 0 2px rgba(0,0,0,.35);pointer-events:none;left:-50px;top:-50px;transform:translate(-50%,-50%)';
(document.body||document.documentElement).appendChild(d);};
addEventListener('mousemove',e=>{mk();const d=document.getElementById('__cur');if(!d)return;d.style.display='block';d.style.left=e.clientX+'px';d.style.top=e.clientY+'px';},true);
document.addEventListener('mouseleave',()=>{const d=document.getElementById('__cur');if(d&&window.top===window)d.style.display='none';});})();"""
def audio(key):
    p=OUT+f"a_{key}.wav"
    if not os.path.exists(p): speak(SEGS[key],p)
    x,sr=sf.read(p); return p,len(x)/sr
class Rec:
    def __init__(s,pw,url,tag):
        s.tag=tag; s.b=pw.chromium.launch(channel="chrome",headless=True)
        s.vd=OUT+"vid_"+tag; shutil.rmtree(s.vd,ignore_errors=True)
        s.ctx=s.b.new_context(viewport={"width":1512,"height":860},record_video_dir=s.vd,record_video_size={"width":1512,"height":860})
        s.ctx.add_init_script(CURSOR); s.pg=s.ctx.new_page(); s.c0=time.time()
        s.pg.goto(BASE+url); s.pg.wait_for_timeout(9000); s.pg.mouse.move(756,430)
        s.t0=time.time()-s.c0; s.segs=[]
    def seg(s,key,fn,lead=0):
        ap,d=audio(key); st=time.time(); s.segs.append((key,ap,st-s.c0-s.t0+lead,d))
        fn(s.pg, lambda: max(0,(st+lead+d+0.5)-time.time()))
        left=(st+lead+d+0.7)-time.time()
        if left>0: s.pg.wait_for_timeout(int(left*1000))
        else: print("  WARN seg overran",key,round(-left,1))
    def finish(s):
        vp=s.pg.video.path(); s.ctx.close(); s.b.close()
        end=s.segs[-1][2]+s.segs[-1][3]+0.7
        cmd=[FF,"-y","-ss",f"{s.t0:.2f}","-t",f"{end:.2f}","-i",vp]
        for _,ap,_,_ in s.segs: cmd+=["-i",ap]
        fl="".join(f"[{i+1}:a]adelay={int(off*1000)}|{int(off*1000)}[a{i}];" for i,(_,_,off,_) in enumerate(s.segs))
        fl+="".join(f"[a{i}]" for i in range(len(s.segs)))+f"amix=inputs={len(s.segs)}:normalize=0[a]"
        cmd+=["-filter_complex",fl,"-map","0:v","-map","[a]","-c:v","libx264","-preset","medium","-crf","20","-pix_fmt","yuv420p","-r","30","-c:a","aac","-b:a","160k",OUT+s.tag+".mp4"]
        r=subprocess.run(cmd,capture_output=True,text=True)
        if r.returncode: print(r.stderr[-1500:]); raise SystemExit(1)
def glide(pg,pts,total,steps=35):
    per=total/len(pts)
    for x,y in pts:
        t=time.time(); pg.mouse.move(x,y,steps=steps); r=per-(time.time()-t)
        if r>0: pg.wait_for_timeout(int(r*1000))
def spin(pg,speed=0.35):
    f=xframe(pg); f.evaluate("(s)=>window.__P.forEach(R=>R.viewer.plugin.canvas3d.setProps({trackball:{animate:{name:s?'spin':'off',params:s?{speed:s}:{}}}}))",speed)
