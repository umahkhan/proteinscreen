import subprocess, os, sys, imageio_ffmpeg, soundfile as sf
from kokoro_onnx import Kokoro
FF=imageio_ffmpeg.get_ffmpeg_exe(); VOICE="af_heart"
V=os.path.expanduser("~/proteinscreen/docs/tutorial/voices/")
K=Kokoro(V+"kokoro-v1.0.onnx",V+"voices-v1.0.bin")
PRON=[("PI3K alpha","P-I-3-K alpha"),("PI3K","P-I-3-K"),("KRAS","K-RAS"),("G12C","G-twelve-C"),("WRN","W-R-N"),("BBO-10203","B-B-O ten two oh three"),("PocketMiner","Pocket Miner"),("ADMET","ad-met"),
      ("GNINA","Gneena"),("9S1A","nine S one A"),("7GQU","seven G Q U"),("9B4U","nine B four U"),("Cys727","cysteine seven twenty-seven"),("Cys242","cysteine two forty-two"),("Cys12","cysteine twelve"),("PyMOL","Pie-mol"),(" A,"," angstroms,"),(" A."," angstroms.")]
def speak(t,path):
    for a,b in PRON: t=t.replace(a,b)
    x,sr=K.create(t,voice=VOICE,speed=0.95,lang="en-us"); sf.write(path,x,sr)
P="pv/"
steps=[
("card_title.png","Covalent drug discovery often starts with a cysteine you can reach, found by chemoproteomics. But the structure you have is usually an apo one, where the pocket may be closed. pscreen is a local pipeline built for that situation. Let's look at what it does on three real targets: WRN helicase, PI3K alpha, and KRAS G12C."),
("wf_preset.png","You start from a preset. This one maps cryptic pockets with two independent models, PocketMiner and Lacuna, names the cysteine you care about, then docks each compound twice: normally, and covalently, bonded to that cysteine. ADMET predictions are optional. The same pipeline also runs from a single command line."),
("runs_wrn.png","Every run keeps a live log with timings. This W-R-N run took about two hours fourteen minutes on a laptop, almost all of it covalent docking, which is slow on a CPU."),
("wrn_overlap.png","First, W-R-N helicase. We start from an apo structure, nine S one A, and cysteine seven twenty-seven, the site targeted by Vividion's covalent W-R-N compounds. The two pocket models map the protein independently. Green is PocketMiner, blue is Lacuna, and purple marks residues where both agree. The yellow marker is the cysteine."),
("wrn_pose.png","Docking Vividion's W-R-N compound covalently, the best pose forms a modeled bond of one point eight angstroms to cysteine seven twenty-seven, inside a pocket that the apo structure does not show open. The legend reports that distance for every compound."),
("pi3k_pose.png","Second, P-I-3-K alpha, a small apo domain, with BridgeBio's B-B-O ten two oh three on cysteine two forty-two. Again the covalent bond is built, at one point eight angstroms. But that is not the whole story."),
("card_check.png","We compared each pose with the real crystal structure. For W-R-N, the top-ranked pose lands about four angstroms from the real ligand, and the closest valid pose about three. For P-I-3-K alpha, the pose sits about eleven angstroms away, so the placement was not reproduced. A modeled bond is not a validated pose, and the tool lets you see both outcomes."),
("kras_pose.png","Third, K-RAS G-twelve-C, a different class of target, from a ligand-free structure, with sotorasib. The covalent pose reaches cysteine twelve, but its score is weak, as you would expect from a closed pocket. Covalent and non-covalent results are ranked in separate lanes and never mixed."),
("card_close.png","So what does this give you? From an apo structure and a cysteine, one reproducible run produces candidate pockets, checked covalent poses, ADMET flags, a report, and a Pie-mol file. What it does not give you is proof. Pockets and poses are hypotheses, scores are weak in closed pockets, covalent docking was limited to sites near the cysteine, and a modeled bond is not evidence of engagement. Use it to decide which experiments to run first."),
]
out=os.path.expanduser("~/proteinscreen/docs/tutorial"); w="pv/_work"; os.makedirs(w,exist_ok=True)
clips=[]
for i,(img,txt) in enumerate(steps):
    wav=f"{w}/a{i}.wav"; speak(txt,wav)
    dur=float(subprocess.check_output(["afinfo",wav]).decode().split("estimated duration:")[1].split()[0])+0.8
    mp4=f"{w}/c{i}.mp4"
    subprocess.run([FF,"-y","-loglevel","error","-loop","1","-i",P+img,"-i",wav,"-t",f"{dur:.2f}","-vf","scale=1512:860,format=yuv420p,fps=10","-c:v","libx264","-c:a","aac","-ar","44100","-ac","2","-af","apad","-shortest",mp4],check=True)
    clips.append((mp4,dur)); print(i,img,round(dur,1),flush=True)
open(w+"/list.txt","w").write("".join(f"file '{os.path.abspath(c[0])}'\n" for c in clips))
subprocess.run([FF,"-y","-loglevel","error","-f","concat","-safe","0","-i",w+"/list.txt","-c","copy",out+"/pscreen_pitch.mp4"],check=True)
print("total",round(sum(c[1] for c in clips)))
