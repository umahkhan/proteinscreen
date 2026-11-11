"""Runs in the isolated PocketMiner env (python3.10 / TF 2.10), NOT the main env.
Usage: pocketminer_infer.py <ae-pocketminer repo> <model: pocketminer|aepocketminer> <in.pdb> <out.csv>
Writes per-residue cryptic-pocket probabilities (chain, resseq, resname, probability).
Code: bowman-lab/ae-pocketminer (MIT); weights: models/pocketminer (original PocketMiner, Meller et al. 2023, MIT)."""
import sys, os
repo, which, pdb, out = sys.argv[1:5]
sys.path.insert(0, os.path.join(repo, "src"))
import numpy as np, tensorflow as tf, mdtraj as md
from models import MQAModel
from util import load_checkpoint
from xtal_predict import process_strucs

att = which == "aepocketminer"
model = MQAModel(node_features=(8, 50), edge_features=(1, 32), hidden_dim=(16, 100),
                 attention_heads=2 if att else None, num_layers=4, dropout=0.1, use_attention=att,
                 save_attention_weights=False, attention_weights_filename=None)
load_checkpoint(model, tf.keras.optimizers.Adam(), os.path.join(repo, "models", which))
s = md.load(pdb)
bb = s.atom_slice(s.top.select("protein and (name N or name CA or name C or name O)"))
X, S, mask = process_strucs([s])
pred = np.asarray(model(X, S, mask, train=False, res_level=True)).reshape(-1)
res = list(bb.top.residues)
assert len(res) == len(pred), (len(res), len(pred))
with open(out, "w") as f:
    f.write("chain,resseq,resname,probability\n")
    for r, p in zip(res, pred):
        f.write(f"{r.chain.chain_id if hasattr(r.chain,'chain_id') and r.chain.chain_id else chr(65 + r.chain.index)},{r.resSeq},{r.name},{float(p):.5f}\n")
