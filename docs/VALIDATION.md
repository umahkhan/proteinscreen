# Retrospective validation — status, results, caveats

Scripts: `scripts_redock.py` / `scripts_bench_all.py` (redocking, covalent), `scripts_cryptic.py` (cryptic pockets); modules in `pscreen/validation/`.
Everything below is retrospective on public data. None of it supports claims about novel targets, potency, selectivity or reaction kinetics.

## Dataset access (verified 2026-10)
| Dataset | Status | Used? |
|---|---|---|
| PDBbind | Site did not respond to scripted access; requires registration under academic terms. Not downloaded. | No (not integrated; supply it yourself after registering) |
| Binding MOAD | Retired in 2024; `bindingmoad.org` now serves an unrelated site. Only an unofficial GitHub mirror of the last affinity CSV exists (licence unclear). | No |
| GalaxyCDock set (`seoklab/cdock_dataset`, `cov_set.json`, 647 covalent complexes) | Public on GitHub; repo has no licence file (PDB entries themselves are public). Data kept local, not redistributed. | Yes (covalent pose recovery) |
| DUD-E | Per-target downloads return HTTP 403 to scripted clients. | No |
| PocketMiner dataset (MIT, bundled in ae-pocketminer) | Public: 38 apo/holo pairs, per-residue labels, val/test split. | Yes (cryptic pockets) |
| RCSB PDB entries | Public. | Yes |

## 1. Redocking (pose recovery), 12 classic noncovalent complexes
Ligand built from CCD SMILES (fresh RDKit conformer), box centred on the crystal ligand, receptor = all chains within 5 Å of the ligand,
exhaustiveness 8, 9 poses, symmetry-aware heavy-atom RMSD, success ≤ 2 Å.
- **v1** (receptor = ligand's own chain only; an input-construction bug that truncated pockets spanning subunits): CNN pose top-1 6/12, Vina-like 4/12.
- **v2** (all chains within 5 Å; fix chosen on principle, not by score): CNN pose top-1 **8/12 (67%)**, Vina-like 6/12 (50%), best-of-9 oracle 10/12 (83%).
  With n = 12 the 95% Wilson interval for 8/12 is roughly 0.39–0.86.
- Failures: sampling (1GPK, 1HQ2: no good pose generated) and scoring (1M17, 1L2S: good pose generated but not ranked first).
- **Leakage caveat:** GNINA's CNN models were trained on PDBbind/CrossDocked-style data; these classic complexes are very likely in it, so results are optimistic.
- Cofactors, waters and ions are removed from receptors (1HQ2's Mg-ATP analogue site is therefore incomplete).

## 2. Covalent pose recovery (GalaxyCDock set), pending completion
Design: ≤30 heavy atoms, one complex per UniProt target (target-level split), 22 selected → 18 buildable (others: mmCIF-only entries or graph mismatch). Ligand = CCD graph with lost leaving atoms removed (form A);
our adduct transform (form B) only where it reproduces the deposited graph. Each ligand docked covalently and noncovalently; success ≤ 2 Å.
Observation already established: GNINA covalent mode returns poses where the bond is **not** satisfied; only poses with attachment atom ≤ 2.3 Å from the nucleophile are valid covalent poses (the app enforces this).
Pose recovery ≠ reaction kinetics, potency or selectivity. **Results: pending (run in progress).**

## 3. Cross-target enrichment (stand-in for DUD-E/PDBbind), pending
Each target's crystal ligand ranked against the other 11 ligands (all docked into every pocket); heavy-atom-count baseline reported beside every signal.
Weak design: decoys are not guaranteed non-binders and not property-matched; one positive per target. **Results: pending.**

## 4. Cryptic-pocket discovery (PocketMiner dataset), TEST split n = 35 (24 with a cryptic pocket, 11 without)
Apo input only; truth = residues labelled 1; hit = Jaccard ≥ 0.25 or centroid ≤ 4 Å; top-k per protein. PocketMiner's residue probability threshold (0.7) was chosen on the validation split and fixed before scoring the test split. Run: `scripts_cryptic.py`, output `runs/cryptic/summary.json`.

| Method | top-1 | top-3 | top-5 | proteins with ≥1 pocket |
|---|---|---|---|---|
| P2Rank 2.5 | 0.524 | 0.714 | 0.762 | 21 |
| Lacuna 1.2.0 | 0.476 | 0.619 | 0.667 | 21 |
| PocketMiner (our clustering of its residue scores) | 0.476 | 0.476 | 0.476 | 21 |

PocketMiner residue-level, pooled over the test set: ROC-AUC 0.868, PR-AUC 0.814 (the paper reports ROC-AUC 0.87, so the local reproduction agrees).
Caveats: (i) the pocket-recovery numbers for PocketMiner depend on **our** residue-clustering step (threshold, 8 Å link, ≥3 residues), which is not part of the published method — they are not a PocketMiner benchmark; (ii) at threshold 0.7, 9 of 11 negative proteins still had at least one residue above threshold, so a bare threshold is a weak "has a cryptic pocket?" classifier; (iii) the three methods answer different questions (residue propensity vs pocket geometry vs ensemble-opened pockets), so the columns are not a ranking of tools; (iv) n = 35 gives wide confidence intervals; (v) Lacuna's claims could only be checked through secondary sources, this table is our own measurement.

## 5. Single-case test: WRN helicase, Vividion covalent inhibitor VVD-133214 (apo input)
Purpose: end-to-end smoke test on a real, recent covalent-allosteric case, not a benchmark. Inputs: apo WRN helicase 517–941 chain A (PDB 9S1A, no nucleotide or DNA); VVD-133214 = RO7589831 (SMILES from Guide to Pharmacology; it matches the 7GQU ligand name); controls: R-epimer, saturated (non-electrophilic) analogue, acrylamide and chloroacetamide fragments, benzamide. Reference: crystal structure of the covalent complex, PDB 7GQU (ligand attached to Cys727 SG, 1.8 Å), superposed on the apo structure by iterative core Cα fitting (core RMSD 0.67 Å; the domain-level shift upon binding is large, 1.7 Å over all common Cα).
See `runs/wrn_apo` (config, manifest, poses). Pocket models: P2Rank, PocketMiner, Lacuna (7 pockets); target residue A:CYS727; 6 compounds.

**Results (single case, no statistics):**
- *Site finding.* Only one of seven predicted pockets is at the Vividion site: P2Rank pocket 3 (its lowest-ranked), centre 2.5 Å from the 7GQU ligand centroid (Cys727 Cα is 5.4 Å from that centroid). The other six are 14.9–24.1 Å away (P2Rank 1 and 2 ≈ 18 Å, PocketMiner 14.9 Å, Lacuna 16.9–24.1 Å). So PocketMiner and Lacuna did not recover this site on this apo structure.
- *Pose recovery (VVD-133214, 30 heavy atoms, MCS-matched RMSD to the crystal pose, apo-frame superposition).* Non-covalent: all five poses 6.9–7.9 Å (top CNN pose score 0.83, Vina-like −5.8 kcal/mol). Covalent (bond to Cys727 SG modelled at 1.81 Å): only 6 of 65 returned poses satisfied the bond; the best of them is 7.4 Å, none within 4 Å. The pose was **not** reproduced.
- *Plausible reasons (untested):* docking into a rigid apo receptor when binding compacts the helicase domain (1.7 Å Cα shift overall); the covalent mode returned many more poses than requested (65 vs `num_modes` 9), which we did not investigate. A holo-receptor (7GQU) re-dock as positive control would separate method limits from apo-structure limits; **not done**.
- *Controls.* Covalent best-pose Vina-like scores: R-epimer −7.0, VVD-133214 −5.2, acrylamide −2.8, chloroacetamide −2.1 kcal/mol; the scores do not reflect the known stereo/warhead preference and must not be read as potency or selectivity. Non-covalent VVD-133214 vs its saturated (non-electrophilic) control at the site: −5.8 vs −5.5, indistinguishable.
- *Run notes.* Covalent docking was restricted to pockets whose centre lies within 12 Å of Cys727 SG (`covalent_max_center_A`) purely for speed; this biases toward the known site and is off by default. The first VVD-133214 covalent job hit the 900 s timeout (the machine was asleep) and was re-run.

## Not performed
Affinity regression (no PDBbind/MOAD access; no comparable labels). Docking into apo vs holo vs Lacuna-open conformations for the cryptic set (planned; not run).
