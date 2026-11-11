from PIL import Image, ImageDraw, ImageFont
W,H=1512,860; LM=100; MAXW=W-2*LM
F="/System/Library/Fonts/Helvetica.ttc"
def font(sz,bold=False): return ImageFont.truetype(F,sz,index=1 if bold else 0)
NAVY,WHITE,MUT,TEAL,AMB="#14143a","#ffffff","#a9adc9","#3fd0c0","#ffb066"
def wrap(d,text,f,maxw):
    out,cur=[],""
    for w in text.split(" "):
        t=(cur+" "+w).strip()
        if d.textlength(t,font=f)<=maxw: cur=t
        else: out.append(cur); cur=w
    return out+[cur]
def card(path,kicker,title,lines,table=None,tsz=54,lsz=32):
    im=Image.new("RGB",(W,H),NAVY); d=ImageDraw.Draw(im); d.rectangle([0,0,14,H],fill=TEAL)
    d.text((LM,100),kicker.upper(),font=font(24,True),fill=TEAL); y=146
    for t in wrap(d,title,font(tsz,True),MAXW): d.text((LM,y),t,font=font(tsz,True),fill=WHITE); y+=int(tsz*1.25)
    y+=22
    for t,col in lines:
        if not t: y+=22; continue
        for s in wrap(d,t,font(lsz),MAXW): d.text((LM,y),s,font=font(lsz),fill=col); y+=int(lsz*1.5)
    if table:
        y+=14; xs=[LM,610,780,1240]
        for i,h in enumerate(table[0]): d.text((xs[i],y),h,font=font(24,True),fill=MUT)
        y+=42; d.line([LM,y-6,W-LM,y-6],fill="#3a3a6a",width=2)
        for row in table[1:]:
            for i,c in enumerate(row):
                col=WHITE if i<3 else (TEAL if row[3]=="close" else AMB)
                d.text((xs[i],y),c,font=font(26,i in (0,3)),fill=col)
            y+=58
    im.save(path)
card("pv/card_title.png","pscreen · pitch","Cryptic pockets + covalent docking on apo structures",
     [("A local, reproducible pipeline for covalent targets where only a closed (apo) structure exists.",WHITE),("",WHITE),
      ("Three real targets:   WRN helicase   ·   PI3K alpha   ·   KRAS G12C",TEAL),("",WHITE),
      ("Public tools and pretrained models only. Predictions to prioritize experiments, not proof of binding.",MUT)],tsz=64)
card("pv/card_check.png","honest check against crystal structures","Is the modeled pose where the real ligand sits?",
     [("Distance between the docked pose centre and the crystal ligand centre (apo frame, backbone-aligned):",MUT)],
     table=[["Target · compound","Crystal ref.","Docked pose (covalent)","Verdict"],
            ["WRN · Cys727 · Vividion compound","7GQU","top 4.2 A · best 3.0 A","close"],
            ["PI3K alpha · Cys242 · BBO-10203","9B4U","10.8 to 11.5 A (6 poses)","not reproduced"],
            ["KRAS G12C · Cys12 · sotorasib","none","no check made","n/a"]])
card("pv/card_close.png","what to take away","From an apo structure and a cysteine to a reproducible shortlist",
     [("In one run: candidate pockets from two independent models, covalent poses with the bond geometry checked, separate covalent and non-covalent rankings, ADMET flags, a report and a PyMOL file.",WHITE),("",WHITE),
      ("Know the limits",AMB),
      ("• Cryptic pockets and docked poses are hypotheses; scores are weak when the pocket is closed.",WHITE),
      ("• Covalent docking was limited to pockets near the target cysteine (17 A for WRN and PI3K alpha, 5 A for KRAS).",WHITE),
      ("• A modeled bond is not evidence of engagement, potency or selectivity.",WHITE),
      ("• The BBO-10203 pose did not match its crystal structure; the VVD-442 structure is reconstructed.",WHITE)],tsz=50,lsz=30)
