# Supplementary Materials

**Author:** Dr Siddalingaiah H S  
**Affiliation:** Professor, Community Medicine, Shridevi Institute of Medical Sciences and Research Hospital, Tumkur  
**Correspondence:** hssling@yahoo.com | +91 8941087719  
**ORCID:** 0000-0002-4771-8285  

## S1. Top pathways
| Pathway | N genes | N significant | Enrichment score |
| :--- | :--- | :--- | :--- |
| mtu00010_Glycolysis | 30 | 29 | 8.6509 |
| mtu00020_TCA_cycle | 25 | 0 | 3.7137 |
| mtu00190_Oxidative_phosphorylation | 35 | 0 | 2.9899 |
| mtu00230_Purine_metabolism | 40 | 0 | 2.8079 |
| mtu00240_Pyrimidine_metabolism | 30 | 0 | 2.4112 |
| mtu00071_Fatty_acid_degradation | 30 | 0 | 2.2621 |
| mtu00640_Propanoate_metabolism | 20 | 0 | 2.2485 |
| mtu00620_Pyruvate_metabolism | 25 | 0 | 2.2328 |
| mtu00550_Peptidoglycan_biosynthesis | 25 | 0 | 2.2145 |
| mtu02020_Two_component_system | 40 | 0 | 1.9935 |

## S2. Top docking results
| Rank | Compound ID | SMILES | Affinity (kcal/mol) | Interacting residues |
| :--- | :--- | :--- | :--- | :--- |
| 1 | MDR_AI_030 | `c1ccccc1-C#C-c1cc(cc(c1)C(F)(F)F)-C#C-C(F)(F)F` | -9.77 | SER-315 |
| 2 | MDR_AI_047 | `c1ccccc1-O-C1C(O)C(O)C(O)C(O)C1-O-C1CCCCC1` | -9.38 | HIS-132, TYR-158 |
| 3 | MDR_AI_026 | `OC-S-c1cc(cnc1)C(=O)N-S-C(F)(F)F` | -8.84 | SER-315, TRP-222, TYR-158 |
| 4 | MDR_AI_010 | `C1CCCCC1-C#C-C1C(O)C(O)C(O)C(O)C1-C#C-C1CCCCC1` | -8.68 | TRP-222 |
| 5 | MDR_AI_039 | `c1ccccc1-O-c1cc(cc(c1)C(F)(F)F)-O-N(C)C` | -8.67 | ASP-94 |
| 6 | MDR_AI_033 | `c1ccccc1-O-c1ccnc2c1cc(cc2)-O-N(C)C` | -8.53 | ASN-426 |
| 7 | MDR_AI_005 | `OC-O-c1cc(cc(c1)C(F)(F)F)-O-C1CCCCC1` | -8.06 | TRP-222 |
| 8 | MDR_AI_019 | `c1ccccc1-O-c1ccnc2c1cc(cc2)-O-c1ccccc1` | -8.03 | ASP-94, TYR-158 |
| 9 | MDR_AI_006 | `c1ccccc1-O-C1C(O)C(O)C(O)C(O)C1-O-N(C)C` | -7.68 | SER-315, TRP-222, TYR-158 |
| 10 | MDR_AI_042 | `OC-C#C-c1cc(cnc1)C(=O)N-C#C-N(C)C` | -7.65 | TYR-158 |

## S3. Model comparison
| Algorithm | Accuracy | Precision | Recall | ROC-AUC | Confusion matrix |
| :--- | :--- | :--- | :--- | :--- | :--- |
| random_forest | 0.650 | 0.600 | 0.667 | 0.727 | [[7, 4], [3, 6]] |
| gradient_boosting | 0.600 | 0.556 | 0.556 | 0.626 | [[7, 4], [4, 5]] |
| logistic_regression | 0.350 | 0.167 | 0.111 | 0.343 | [[6, 5], [8, 1]] |

## S4. Mutation catalog excerpt
| Gene | Mutation | Drug | Frequency | Confidence |
| :--- | :--- | :--- | :--- | :--- |
| Rv0678 | V1F | Bedaquiline/Clofazimine | 0.10 | moderate |
| Rv0678 | M1R | Bedaquiline/Clofazimine | 0.08 | moderate |
| Rv0678 | R109Stop | Bedaquiline/Clofazimine | 0.05 | moderate |
| atpE | E61D | Bedaquiline | 0.15 | high |
| atpE | A63P | Bedaquiline | 0.10 | high |
| atpE | I66M | Bedaquiline | 0.08 | high |
| atpE | D28V | Bedaquiline | 0.05 | high |
| embB | M306V | Ethambutol | 0.35 | moderate |
| embB | M306I | Ethambutol | 0.15 | moderate |
| embB | G406A | Ethambutol | 0.08 | moderate |
| embB | Q497R | Ethambutol | 0.06 | moderate |
| ethA | A381P | Ethionamide | 0.12 | moderate |
| ethA | Q254Stop | Ethionamide | 0.08 | moderate |
| ethA | Y84D | Ethionamide | 0.05 | moderate |
| gyrA | D94G | Fluoroquinolones | 0.30 | high |
| gyrA | A90V | Fluoroquinolones | 0.20 | high |
| gyrA | D94A | Fluoroquinolones | 0.10 | high |
| gyrA | D94N | Fluoroquinolones | 0.08 | high |
| gyrA | D94Y | Fluoroquinolones | 0.05 | high |
| inhA | C-15T | Isoniazid (low-level) | 0.25 | high |
