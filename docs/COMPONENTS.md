# Component status and limitations (updated as work proceeds)

| Component | Status | License (code / weights) | Notes |
|---|---|---|---|
| RDKit prep, dedup, descriptors, warhead SMARTS (17 classes) | Implemented, tested | BSD | Warhead match = annotation only. Reactive-atom indices refer to the prepared ligand SDF heavy-atom order, not canonical SMILES. |
| GNINA noncovalent docking (CPU, via Lima+Rosetta) | Implemented, smoke-tested | Apache-2.0/GPL | ~20-40 s per small ligand (3 CPUs). Large boxes/ligands slower. |
| GNINA covalent docking | NOT yet implemented/verified | | Flags exist (`--covalent_*`); supported chemistries untested. |
| Reactive-atom/target geometry (noncovalent poses) | Implemented | | Distance thresholds (4/6 A) are our assumptions. No attack-angle check yet. |
| P2Rank | Implemented | MIT | Conventional pockets only. |
| PocketMiner | Implemented | MIT / MIT (checkpoint bundled in ae-pocketminer) | Per-residue probabilities only. Threshold 0.7, 8 A single-linkage grouping, ranking by summed probability are OUR choices, unvalidated; on dimers/compact sites it can merge residues into one large cluster (box capped at 32 A). Checkpoint provenance relative to the original repo not byte-compared. |
| Lacuna | Implemented | MIT | NMA ensemble is a coarse model; its `cryptic` flag also marks known open sites. Benchmarks are the authors' own preprint. Docking targets the largest-volume conformer. |
| Chemprop | Not started | MIT | No general pretrained checkpoints ship with the package. |
| Validation harness | Not started | | Needs symmetry-aware RMSD (spyrmsd/RDKit). |
