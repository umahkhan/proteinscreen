# End-to-end test results

Each row is one run of `python -m pscreen run --fast ...` (GNINA exhaustiveness 4, 3 poses, CPU, 6-core VM). Scores are docking predictions, not affinities.

| run          | result   |   seconds | sites                                                                |   compounds in |   docked |   not docked | top compound per lane                                                                               |
|:-------------|:---------|----------:|:---------------------------------------------------------------------|---------------:|---------:|-------------:|:----------------------------------------------------------------------------------------------------|
| brd4         | SUCCESS  |       142 | p2rank: ok (2 pockets); pocketminer: ran, no pockets above threshold |              4 |        3 |            1 | noncovalent: JQ1 (VS 3.86)                                                                          |
| caii         | SUCCESS  |        41 | user_box: ok                                                         |              3 |        3 |            0 | noncovalent: benzenesulfonamide (VS 6.11)                                                           |
| hivpr        | SUCCESS  |        85 | known_ligand: ok                                                     |              3 |        3 |            0 | noncovalent: indinavir (VS 5.54)                                                                    |
| kras         | SUCCESS  |       621 | p2rank: ok (3 pockets)                                               |              3 |        3 |            0 | covalent: acryloyl_phenylpiperazine (VS 2.45); noncovalent: chloroacetyl_phenylpiperazine (VS 3.02) |
| streptavidin | SUCCESS  |        67 | p2rank: ok (1 pockets)                                               |              4 |        3 |            1 | noncovalent: biotin (VS 8.90)                                                                       |
| t4l          | SUCCESS  |       219 | lacuna: ok (3 pockets)                                               |              4 |        4 |            0 | noncovalent: toluene (VS 3.89)                                                                      |
| trypsin      | SUCCESS  |        54 | user_box: ok                                                         |              3 |        3 |            0 | noncovalent: 4-aminobenzamidine (VS 3.91)                                                           |

## Per-compound results

| run          | lane        | rank   | compound                      | site         | CNN VS   | kcal/mol   | warhead         | cov. bond Å                                            |
|:-------------|:------------|:-------|:------------------------------|:-------------|:---------|:-----------|:----------------|:-------------------------------------------------------|
| brd4         | noncovalent | 1      | JQ1                           | p2rank_1     | 3.86     | -8.99      |                 |                                                        |
| brd4         | noncovalent | 2      | glucose_decoy                 | p2rank_2     | 2.99     | -3.63      |                 |                                                        |
| brd4         | noncovalent | 3      | dimethylisoxazole             | p2rank_1     | 2.72     | -3.49      |                 |                                                        |
| brd4         | —           |        | JQ1_duplicate                 |              |          |            |                 | duplicate: duplicate of JQ1                            |
| caii         | noncovalent | 1      | benzenesulfonamide            | user_box     | 6.11     | -6.46      |                 |                                                        |
| caii         | noncovalent | 2      | acetazolamide                 | user_box     | 5.7      | -6.17      |                 |                                                        |
| caii         | noncovalent | 3      | ethanol_decoy                 | user_box     | 2.78     | -2.04      |                 |                                                        |
| hivpr        | noncovalent | 1      | indinavir                     | known_ligand | 5.54     | -10.02     |                 |                                                        |
| hivpr        | noncovalent | 2      | amprenavir_like_fragment      | known_ligand | 5.33     | -8.32      |                 |                                                        |
| hivpr        | noncovalent | 3      | caffeine_decoy                | known_ligand | 3.38     | -5.23      |                 |                                                        |
| kras         | covalent    | 1      | acryloyl_phenylpiperazine     | p2rank_1     | 2.45     | -3.09      | acrylamide      | 1.81                                                   |
| kras         | covalent    | 2      | chloroacetyl_phenylpiperazine | p2rank_1     | 2.09     | -2.9       | chloroacetamide | 1.81                                                   |
| kras         | noncovalent | 1      | chloroacetyl_phenylpiperazine | p2rank_2     | 3.02     | -3.05      | chloroacetamide | 7.56                                                   |
| kras         | noncovalent | 2      | propanoyl_control             | p2rank_2     | 2.63     | -3.09      |                 |                                                        |
| kras         | noncovalent | 3      | acryloyl_phenylpiperazine     | p2rank_1     | 2.35     | -3.06      | acrylamide      | 3.92                                                   |
| streptavidin | noncovalent | 1      | biotin                        | p2rank_1     | 8.9      | -7.55      |                 |                                                        |
| streptavidin | noncovalent | 2      | desthiobiotin                 | p2rank_1     | 6.8      | -7.5       |                 |                                                        |
| streptavidin | noncovalent | 3      | phenol_decoy                  | p2rank_1     | 3.39     | -2.73      |                 |                                                        |
| streptavidin | —           |        | not_a_smiles                  |              |          |            |                 | invalid: SMILES could not be parsed/sanitized by RDKit |
| t4l          | noncovalent | 1      | toluene                       | lacuna_2     | 3.89     | -6.24      |                 |                                                        |
| t4l          | noncovalent | 2      | benzene                       | lacuna_2     | 3.57     | -5.46      |                 |                                                        |
| t4l          | noncovalent | 3      | glycine_decoy                 | lacuna_1     | 3.31     | -3.39      |                 |                                                        |
| t4l          | noncovalent | 4      | indole                        | lacuna_2     | 3.13     | -6.15      |                 |                                                        |
| trypsin      | noncovalent | 1      | 4-aminobenzamidine            | user_box     | 3.91     | -6.21      |                 |                                                        |
| trypsin      | noncovalent | 2      | benzamidine                   | user_box     | 3.62     | -6.06      |                 |                                                        |
| trypsin      | noncovalent | 3      | benzene_decoy                 | user_box     | 1.72     | -3.37      |                 |                                                        |
