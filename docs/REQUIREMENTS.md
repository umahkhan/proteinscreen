---
title: "pscreen: Requirements Document"
---

**Version** 0.1 · **Status** Working prototype, public pilot hosted on a personal computer · **Document type** Product requirements

# 1. Purpose

pscreen is a local, reproducible tool that helps scientists prioritize which compounds to test first against a protein target. It combines binding-pocket discovery (including cryptic pockets that are closed in the available structure), compound docking (including covalent docking to a chosen cysteine), compound annotation, and interactive 3D review in one workflow that needs no programming.

# 2. Problem statement

Structure-based screening currently means chaining several research tools with different environments, input formats and outputs. Three problems follow:

1. **Access.** Bench scientists without command-line experience cannot run the tools, so computational work gets delegated or skipped.
2. **Cryptic pockets.** Many targets are only available as closed (apo) structures. Standard pocket finders and docking tools miss pockets that open dynamically, so covalent drug discovery on such targets is hard to explore.
3. **Reproducibility and honesty.** Scripted pipelines often lose track of tool versions, settings and failures, and docking scores are easily over-read as evidence of binding.

# 3. Goals and non-goals

**Goals**

- G1. Let a non-programmer go from a protein structure and a compound list to a ranked, explorable shortlist in one session.
- G2. Support apo structures through independent cryptic-pocket predictors, shown side by side and never merged into one score.
- G3. Support covalent docking to a user-specified residue, with the bond geometry checked.
- G4. Make every run reproducible and every failure visible.
- G5. State the limits of the predictions inside the product itself.

**Non-goals**

- Training or fine-tuning models. Only public tools and pretrained models are used.
- Claims about binding affinity, potency, selectivity or covalent engagement. Docking and pocket outputs are hypotheses for prioritizing experiments.
- Replacing expert review or experimental validation.
- Screening libraries beyond about 1,000 compounds in the local version.

# 4. Users

| User | Need | Skill level |
|:---|:---|:---|
| Bench or medicinal chemist | Rank candidate compounds for a known target and see where they sit | Low computational |
| Structural or computational biologist | Explore pockets and compare predictors on an apo structure | High |
| Evaluator or collaborator | Open a result and understand it without installing anything | Low |

# 5. Scope and user flow

1. **Workflow.** Build a pipeline on a node canvas (or pick a preset): protein, compounds, binding-site models, covalent target residue, docking, properties, results.
2. **Runs.** Start the experiment and watch per-step progress and timers.
3. **Results.** Review a ranked table and the 3D explorer, then export.

# 6. Functional requirements

Priority: **M** must, **S** should, **C** could.

## Inputs

- FR-1 (M) Accept a protein as PDB or mmCIF. Remove waters, ions and ligands for docking.
- FR-2 (M) Accept compounds in one standard format: a CSV with a header row `id,smiles`. Also accept TSV, Excel, `.smi` and SDF, and convert them to the standard form.
- FR-3 (M) Provide a downloadable compound template and a preview showing how every row was read. Rows that cannot be used are listed with a reason and are never silently dropped. Examples: invalid SMILES, duplicate id.
- FR-4 (S) Accept an optional known ligand or box to define a site, and optional experimental data to join to results.

## Pocket discovery

- FR-5 (M) Run conventional pocket prediction (P2Rank).
- FR-6 (M) Run cryptic-pocket predictors: a conformational ensemble method (Lacuna) and per-residue cryptic propensity (PocketMiner). Show each tool's own scores. Do not merge scores across tools.
- FR-7 (S) Highlight residues claimed by more than one tool.

## Docking and annotation

- FR-8 (M) Dock compounds with GNINA and report its scores separately (CNN VS, CNN pose, CNN affinity, Vina-like). Never combine them.
- FR-9 (M) Covalent lane: for a user-specified residue, model the adduct for recognized warheads and report the bond geometry. Keep covalent and non-covalent rankings separate.
- FR-10 (S) Annotate electrophilic warheads (RDKit) and, optionally, ADMET properties (ADMET-AI).

## Workflow and execution

- FR-11 (M) Every combination of selected steps must run as built. Presets only place nodes on the canvas.
- FR-12 (M) Runs are resumable. Finished docking jobs are skipped on restart. Failure of an optional step is recorded and does not stop the run.
- FR-13 (M) Record a manifest per run: tool versions, input hashes and configuration.

## Review and export

- FR-14 (M) A 3D explorer showing pockets (colored by tool), docked poses, the target residue and the warhead-to-residue distance.
- FR-15 (M) A ranked table per docking lane, with selectable ranking signal and a per-compound detail view.
- FR-16 (M) Export a package: report, ranked CSV, per-compound status CSV, pockets CSV, best poses SDF and a PyMOL scene.
- FR-17 (S) A pockets-only mode (no docking) for fast structure triage.

# 7. Scientific-integrity requirements

- SI-1 No silent method substitution. If a tool fails, show the failure and do not replace it with another tool.
- SI-2 Every output view and export carries the statement that results are predictions, not evidence of binding, engagement, potency or selectivity.
- SI-3 Document known limitations in the product. Examples: weak docking scores in closed pockets; a covalent-docking distance limit that biases results toward the target site; reconstructed reference structures.
- SI-4 Cite the underlying tools in the product and exports.

# 8. Non-functional requirements

| ID | Requirement |
|:---|:---|
| NFR-1 Usability | A first-time user completes a run with only the on-screen guidance and the compound template. |
| NFR-2 Reproducibility | The same inputs, settings and seed reproduce the same results on the same tool versions. |
| NFR-3 Performance | Pocket-only runs finish in seconds to minutes. A small docking run (3–4 compounds, two sites) finishes in about 3–4 minutes on a laptop. A full covalent run with all three pocket models took over 2 hours in testing. |
| NFR-4 Privacy | Inputs stay on the machine running the app. The local version sends nothing off the machine. |
| NFR-5 Licensing | Use only tools whose licenses permit the intended use. Do not redistribute third-party binaries or model weights in the repository. |

# 9. Hosted-use requirements (public pilot)

- HR-1 Each visitor gets a private workspace addressed by an unguessable token. Visitors cannot see each other's runs.
- HR-2 Enforce limits on the server side: at most 25 compounds, 15 MB protein, 2 covalent residues, one active run per visitor, six runs queued server-wide. Queue runs on one slot at a time.
- HR-3 Delete workspaces after 48 hours idle.
- HR-4 Sanitize uploaded file names so they cannot write outside the visitor's folder.
- HR-5 Tell visitors that the server operator can read uploaded files, and that confidential material should not be uploaded unless the operator is trusted.
- HR-6 Ship a pockets-only variant suitable for a free managed host with limited memory.

# 10. Constraints and assumptions

- Public tools and pretrained models only. No training.
- GNINA ships as a Linux x86-64 binary. On macOS it runs in a virtual machine.
- PocketMiner needs an older Python and TensorFlow, so it runs in its own isolated environment.
- Free managed hosts do not provide enough compute for docking. The public pilot therefore runs on a personal computer through a temporary tunnel.

# 11. Risks and mitigations

| Risk | Mitigation |
|:---|:---|
| Predictions are over-read as proof | Built-in disclaimers (SI-2), separate score columns, documented limits (SI-3) |
| Uploaded files exploit a parser on the host | Per-visitor isolation, file-name sanitizing, caps; advise a dedicated account or machine for wider audiences |
| Host unavailable (sleep, network) | Document that the pilot is best-effort; plan for paid hosting if demand justifies it |
| Misconfigured or failing tool | Clear per-tool error messages and a manifest; no fallback to a different method |
| License terms change or are unclear | Review third-party licenses before each release; keep a third-party notices file |

# 12. Acceptance criteria and validation

- AC-1 A non-programmer completes a run from the template and the on-screen guidance alone.
- AC-2 On a known system (BRD4 with its reference binder against decoys), the known binder ranks first. *Observed in testing.*
- AC-3 On the WRN helicase apo structure, pocket predictors and the covalent lane report sites and bond geometry near the target cysteine. *Observed in testing.*
- AC-4 A failing step shows a plain-language error, and the rest of the run continues where possible.
- AC-5 In hosted mode, one visitor cannot read another's workspace, and a hostile file name cannot escape its folder. *Observed in testing.*
- AC-6 An automated test suite passes (89 tests at the time of writing).

# 13. Success measures (targets, not yet measured)

- Share of first-time users who complete a run unaided.
- Median time from upload to ranked shortlist.
- Share of runs that finish without a hard failure.
- Number of distinct users and repeat users during the pilot.

# 14. Open questions

1. Which hosting model should follow the pilot: a paid managed host, an institutional server, or local installs only?
2. Which license applies to the project, and do all third-party tools permit hosted public use? One tool's license file needs confirmation.
3. How large should the free-tier limits be, given the cost of docking compute?
4. How should the product support users who need confidentiality: local-only guidance, or an isolated hosted option?
5. What external validation (blinded comparison against known structures) is needed before wider promotion?

# 15. Release status

Local application, command-line interface, hosted mode with private workspaces and a queue, and a pockets-only app are built and tested on a single development machine. The public pilot, the container image for managed hosting and the free managed-host deployment have not been fully verified.
