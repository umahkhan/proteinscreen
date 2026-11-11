# pscreen — pocket discovery and compound prioritization

Upload a protein and up to ~1000 compounds, discover or define binding pockets, dock with GNINA, annotate electrophilic warheads,
review ranked compounds and inspect the real docked poses (3Dmol, Mol*, or PyMOL renders), export everything. See `docs/SETUP.md`
(install), `docs/DEPLOY.md` and `docs/SERVE_FROM_MAC.md` (sharing it with others), `docs/COMPONENTS.md` (what is implemented / limited) and `docs/VALIDATION.md` (benchmarks and their caveats).

## Demo

[![pscreen demo — click to watch](https://img.shields.io/badge/▶_Watch-the_demo_video-blue?style=for-the-badge)](docs/tutorial/pscreen_demo.mp4)

The full walkthrough is at [`docs/tutorial/pscreen_demo.mp4`](docs/tutorial/pscreen_demo.mp4) (opens in GitHub's video player).

Docking scores are computational predictions. They do not establish binding, covalent engagement, potency or selectivity.
A warhead match is a chemical annotation, not proof of reactivity.

## Quick start (CLI — same pipeline as the web app)

```bash
python -m pscreen check --protein 3MXF.pdb --compounds jq1.csv                       # validate inputs only (fast, clear errors)
python -m pscreen run   --name brd4 --protein 3MXF.pdb --compounds jq1.csv --pockets p2rank,pocketminer --admet --fast
python -m pscreen run   --name kras --protein 4LDJ_A.pdb --compounds acrylamides.csv --residue A:CYS12 --pockets p2rank   # + covalent docking at Cys12
python -m pscreen run   ... --known-ligand ligand.pdb          # site from a crystal ligand      |   --box=X,Y,Z,SX,SY,SZ  # manual box
python -m pscreen list | export runs/<name>                    # list runs / (re)write the output package
```

**Inputs**: protein (PDB/mmCIF) · compounds (`.csv/.tsv` with a SMILES column, `.smi`, `.sdf`) · optional covalent target residue(s) · optional known ligand or box · optional pocket predictors (`p2rank`, `pocketminer`, `lacuna`) · optional ADMET-AI.
**Outputs** (`runs/<name>/outputs/`): `REPORT.md` (human summary), `ranked_compounds.csv` (one row per compound and docking lane), `compound_status.csv` (every input compound: docked / invalid / duplicate / failed + reason),
`pockets.csv`, `best_poses.sdf`, `summary.json`. The Results page offers the same package as one zip (Export tab). `--fast` = exhaustiveness 4, 3 poses. Test results: `docs/E2E_RESULTS.md` (`python scripts/e2e.py`).

**Screens**: *Workflow* — a node canvas. Left: the workflow steps in pipeline order (inputs → binding-site models → covalent target residue(s) → docking → compound properties → results); expand a step's ▾ to set its parameters, then drag it (or press "Add to canvas"). Middle: the canvas. Right: inputs (experiment name, protein, compound library incl. Excel/CSV/SDF, optional known ligand and experimental data) and **Run analysis**. Any combination is valid, e.g. P2Rank only with no ADMET, or P2Rank + PocketMiner + Lacuna with covalent docking on several residues and ADMET-AI. Below: all experiments (status, current step, elapsed) and a per-step log with a timer for every step plus the total. Each run is its own process, so several workflows can run at the same time (they share the machine's CPUs). *Results* — the ranked-compound table and Mol* 3D explorer for any finished run.