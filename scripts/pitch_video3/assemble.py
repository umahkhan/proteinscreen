import subprocess, re, imageio_ffmpeg
FF=imageio_ffmpeg.get_ffmpeg_exe(); names=["g0a","g0b","g1","g2","g3"]; X=0.6
def dur(f):
    o=subprocess.run([FF,"-i",f],capture_output=True,text=True).stderr; h,m,s=re.search(r"Duration: (\d+):(\d+):([\d.]+)",o).groups(); return int(h)*3600+int(m)*60+float(s)
D=[dur(f"pv/{n}.mp4") for n in names]; cmd=[FF,"-y","-loglevel","error"]
for n in names: cmd+=["-i",f"pv/{n}.mp4"]
fv="";fa="";acc=D[0];pv="[0:v]";pa="[0:a]"
for i in range(1,len(names)):
    fv+=f"{pv}[{i}:v]xfade=transition=fade:duration={X}:offset={acc-X:.2f}[v{i}];"; pv=f"[v{i}]"
    fa+=f"{pa}[{i}:a]acrossfade=d={X}[a{i}];"; pa=f"[a{i}]"; acc+=D[i]-X
cmd+=["-filter_complex",(fv+fa).rstrip(";"),"-map",pv,"-map",pa,"-c:v","libx264","-crf","20","-preset","medium","-pix_fmt","yuv420p","-c:a","aac","-b:a","160k","../../docs/tutorial/pscreen_pitch_v4.mp4"]
subprocess.run(cmd,check=True); print("durations",[round(d,1) for d in D],"total",round(acc,1))
