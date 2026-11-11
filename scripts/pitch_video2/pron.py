"""Pronunciation layer for Kokoro: known terms get explicit IPA; everything else uses espeak. Text -> phoneme string -> is_phonemes=True."""
import re, numpy as np, soundfile as sf
from kokoro_onnx import Kokoro
from kokoro_onnx.tokenizer import Tokenizer
V="/Users/umairkhan/proteinscreen/docs/tutorial/voices/"
K=Kokoro(V+"kokoro-v1.0.onnx",V+"voices-v1.0.bin"); TOK=Tokenizer(); VOICE="af_heart"
# term -> IPA (stress marks: ˈ primary)
LEX={
 "apo":"ˈeɪpoʊ", "cryptic":"kɹˈɪptɪk", "covalent":"koʊvˈeɪlənt", "covalently":"koʊvˈeɪləntli", "cysteine":"sˈɪstiːn",
 "helicase":"hˈɛlɪkeɪs", "Vividion":"vɪvˈɪdiən", "WRN":"dˈʌbəljuː ˈɑːɹ ˈɛn", "Lacuna":"lɐkˈuːnə",
 "PocketMiner":"pˈɑːkɪt mˈaɪnɚ", "GNINA":"nˈiːnə", "ADMET":"ˈæd mˈɛt", "PyMOL":"pˈaɪmoʊl", "RDKit":"ˈɑːɹ dˈiː kˈɪt",
 "Mol*":"mˈoʊl stˈɑːɹ", "angstroms":"ˈæŋstɹəmz", "ligand":"lˈɪɡənd", "ligands":"lˈɪɡəndz", "chemoproteomics":"kˌiːmoʊpɹˌoʊtiːˈɑːmɪks",
 "druggable":"dɹˈʌɡəbəl", "conformer":"kənfˈɔːɹmɚ", "conformers":"kənfˈɔːɹmɚz", "electrophile":"ɪlˈɛktɹəfaɪl", "warhead":"wˈɔːɹhɛd",
 "pscreen":"pˈiː skɹˈiːn", "CLI":"sˈiː ˈɛl ˈaɪ", "SMILES":"smˈaɪlz", "PDB":"pˈiː dˈiː bˈiː", "Cys727":"sˈɪstiːn sˈɛvən twˈɛntisˈɛvən",
}
_re=re.compile(r"(?<![\w*])("+"|".join(sorted(map(re.escape,LEX),key=len,reverse=True))+r")(?![\w])")
def ph(text):
    out=[]; pos=0
    for m in _re.finditer(text):
        seg=text[pos:m.start()]
        if seg.strip(): out.append(TOK.phonemize(seg,"en-us"))
        out.append(LEX[m.group(1)]); pos=m.end()
    if text[pos:].strip(): out.append(TOK.phonemize(text[pos:],"en-us"))
    return " ".join(out)
def speak(text,path=None,speed=0.95):
    x,sr=K.create(ph(text),voice=VOICE,speed=speed,is_phonemes=True)
    if path: sf.write(path,x,sr)
    return x,sr
