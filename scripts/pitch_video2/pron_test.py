from pron import *
words=["apo","cryptic","covalent","covalently","cysteine","helicase","Vividion","WRN","Lacuna","PocketMiner","GNINA","ADMET","PyMOL","RDKit","Mol*","angstroms","ligand","chemoproteomics","druggable","conformer","electrophile","warhead","pscreen","SMILES","PDB"]
sr=24000; sil=np.zeros(int(.55*sr),dtype="float32"); parts=[]
for w in words:
    print(w, "->", ph(w))
    a,sr=speak(w+".",speed=0.9); parts += [a,sil]
a,_=speak("From an apo structure, and cysteine seven twenty-seven on WRN helicase, pscreen finds cryptic pockets and docks covalent ligands.")
parts.append(a)
sf.write("pron_test.wav",np.concatenate(parts),sr)
