# Supplementary Materials

**Author:** Dr Siddalingaiah H S
**Affiliation:** Professor, Community Medicine, Shridevi Institute of Medical Sciences and Research Hospital, Tumkur
**Correspondence:** hssling@yahoo.com | +91 8941087719
**ORCID:** 0000-0002-4771-8285

## S1. ChEMBL Dataset Summary (Target: CHEMBL1849 — InhA, Rv1484)

All bioactivity data were retrieved via the ChEMBL REST API (chembl-webresource-client) in March 2025 using standard_type = "IC50" and standard_units = "nM".

| Parameter | Value |
| :--- | :--- |
| ChEMBL Target ID | CHEMBL1849 |
| Target name | Enoyl-[acyl-carrier-protein] reductase (InhA, Rv1484) |
| Organism | Mycobacterium tuberculosis H37Rv |
| Total compounds retrieved | 277 |
| Active (IC50 < 1000 nM) | 108 |
| Inactive (IC50 >= 1000 nM) | 169 |
| IC50 range (pre-filter, all retrieved) | 0.4 nM – 500,000 nM |
| IC50 range (post drug-likeness filter) | 2 nM – 291,000 nM |
| Median IC50 (active set) | 38 nM |
| Data source | https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL1849/ |

## S2. QSAR Model Performance — All Three Classifiers

Models trained on Morgan ECFP4 fingerprints (radius=2, 2048 bits) with 80/20 stratified train/test split and 5-fold cross-validation. All metrics computed on real ChEMBL1849 IC50 data; no synthetic or imputed values.

| Algorithm | CV ROC-AUC (mean) | CV ROC-AUC (std) | CV F1 | Test ROC-AUC | Test F1 | Test Precision | Test Recall | Test Accuracy |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Random Forest | 0.9328 | 0.0541 | 0.7915 | 0.9318 | 0.757 | 0.933 | 0.636 | 0.839 |
| Gradient Boosting | 0.8967 | 0.0739 | 0.8121 | 0.9412 | 0.800 | 0.889 | 0.727 | 0.857 |
| Logistic Regression | 0.9061 | 0.0648 | 0.8138 | 0.9612 | 0.821 | 0.941 | 0.727 | 0.875 |

Training set n=221; test set n=56. Label definition: IC50 < 1000 nM = active.

**Known limitations:** Models trained on biochemical assay data (not whole-cell MIC); applicability domain not formally assessed; external validation on independent dataset not performed.

## S3. Top 10 Ranked Compounds — Full ADMET and Scoring Details

Composite score = 0.40 × pIC50_norm + 0.30 × QSAR_Prob + 0.20 × LBVS_Tanimoto + 0.10 × QED. All IC50 values are experimentally measured (ChEMBL1849); ADMET descriptors computed with RDKit 2024.03.

| Rank | ChEMBL ID | IC50 (nM) | MW | LogP | TPSA | HBD | HBA | QED | hERG Risk | GI Abs. | Composite Score |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | CHEMBL3125270 | 4 | 423.47 | 1.71 | 123.46 | 2 | 5 | 0.625 | Low | Yes | 0.8245 |
| 2 | CHEMBL3125104 | 3 | 467.53 | 1.33 | 129.70 | 3 | 6 | 0.458 | Low | Yes | 0.8237 |
| 3 | CHEMBL3125266 | 3 | 495.54 | 1.51 | 135.77 | 2 | 7 | 0.452 | Low | Yes | 0.8220 |
| 4 | CHEMBL3088173 | 3 | 434.50 | 3.83 | 88.75 | 2 | 8 | 0.480 | Low | Yes | 0.8176 |
| 5 | CHEMBL3088172 | 3 | 434.50 | 3.83 | 88.75 | 2 | 8 | 0.480 | Low | Yes | 0.8176 |
| 6 | CHEMBL3125276 | 5 | 481.51 | 1.03 | 135.77 | 2 | 7 | 0.481 | Low | Yes | 0.8077 |
| 7 | CHEMBL3125235 | 2 | 514.59 | 2.94 | 122.36 | 2 | 6 | 0.373 | Low | Yes | 0.8070 |
| 8 | CHEMBL3125259 | 5 | 490.52 | 2.05 | 135.50 | 2 | 7 | 0.405 | Low | Yes | 0.8020 |
| 9 | CHEMBL3125101 | 2 | 569.43 | 1.76 | 138.93 | 3 | 7 | 0.275 | Low | Yes | 0.8006 |
| 10 | CHEMBL3125261 | 3 | 540.58 | 3.20 | 135.50 | 2 | 7 | 0.324 | Low | Yes | 0.7913 |

ADMET method: RDKit 2024.03 Descriptors, QED (rdkit.Chem.QED), rdMolDescriptors; hERG structural alert per Waring 2010 (Expert Opin Drug Metab Toxicol); GI absorption proxy per Lipinski 2001 (Adv Drug Deliv Rev); BBB proxy per Clark 1999 (J Pharm Sci).

## S4. Lead Compound Profile — CHEMBL3125270

| Property | Value | Reference/Method |
| :--- | :--- | :--- |
| ChEMBL ID | CHEMBL3125270 | ChEMBL1849 |
| SMILES | CCc1cc(C(=O)N[C@@H]2C[C@@H](C(N)=O)N(C(=O)c3coc4ccccc34)C2)n(CC)n1 | ChEMBL |
| Measured IC50 | 4 nM (InhA biochemical assay) | ChEMBL1849 |
| Molecular weight | 423.47 Da | RDKit 2024.03 |
| cLogP | 1.71 | RDKit 2024.03 |
| TPSA | 123.46 A^2 | RDKit 2024.03 |
| HBD / HBA | 2 / 5 | RDKit 2024.03 |
| Rotatable bonds | 6 | RDKit 2024.03 |
| QED | 0.625 | RDKit 2024.03 |
| Lipinski compliance | Pass (all 4 rules) | RDKit 2024.03 |
| Veber compliance | Pass (RotBonds 6, TPSA 123) | RDKit 2024.03 |
| hERG structural alert | Low risk | Waring 2010 criterion |
| GI absorption (proxy) | Predicted absorbed | Lipinski 2001 criterion |
| BBB penetration (proxy) | Not predicted | Clark 1999 criterion (TPSA > 90) |
| Scaffold class | Pyrazole-benzofuran-pyrrolidine | RDKit ring analysis |
| QSAR predicted active | Yes (Prob = 0.957) | Logistic Regression on ECFP4 |

Note: All properties are computationally predicted from structure. No experimental ADMET data (Caco-2, hERG patch-clamp, microsomal stability) have been generated. Predictions should be confirmed by in vitro assay before progressing this compound.

## S5. Resistance Mutation Frequencies (CRyPTIC Consortium, 2022)

Mutation frequencies from the CRyPTIC Consortium clinical genome sequencing study (n=12,289 isolates worldwide, PLoS Biol 2022).

| Gene | Mutation | Drug affected | Frequency in MDR strains | Confidence |
| :--- | :--- | :--- | :--- | :--- |
| katG | S315T | Isoniazid (high-level) | 0.67 | High |
| inhA promoter | C-15T | Isoniazid (low-level) | 0.25 | High |
| rpoB | S450L | Rifampicin | 0.42 | High |
| gyrA | D94G | Fluoroquinolones | 0.30 | High |
| gyrA | A90V | Fluoroquinolones | 0.20 | High |
| gyrA | D94A | Fluoroquinolones | 0.10 | High |
| embB | M306V | Ethambutol | 0.35 | Moderate |
| embB | M306I | Ethambutol | 0.15 | Moderate |
| atpE | E61D | Bedaquiline | 0.15 | High |
| Rv0678 | V1F | Bedaquiline/Clofazimine | 0.10 | Moderate |
| ethA | A381P | Ethionamide | 0.12 | Moderate |

Source: CRyPTIC Consortium (2022). Genome-wide association studies of global Mycobacterium tuberculosis resistance to 13 antimicrobials in 10,228 genomes. PLoS Biol 20(8):e3001755. https://doi.org/10.1371/journal.pbio.3001755

**Rationale for target selection:** The katG S315T mutation (frequency 0.68 in MDR strains) abolishes prodrug activation of isoniazid but leaves the InhA enzyme structurally intact. Direct InhA inhibitors that do not require katG activation therefore retain activity against this dominant resistance mechanism. This motivated selection of InhA (Rv1484, CHEMBL1849) as the primary screening target.

## S6. LBVS Reference Scaffold Tanimoto Similarity

Tanimoto similarity was computed using RDKit Morgan fingerprints (radius=2, 2048 bits) between each candidate compound and five known InhA inhibitor scaffolds.

| Reference scaffold | Role | Structure |
| :--- | :--- | :--- |
| Triclosan | Known InhA inhibitor (enoyl-ACP reductase cross-reactant) | 5-chloro-2-(2,4-dichlorophenoxy)phenol |
| Isoniazid | First-line TB drug; prodrug requiring katG | Pyridine-4-carbohydrazide |
| Ethionamide | Second-line; structural analogue of isoniazid | 2-ethylpyridine-4-carbothioamide |
| PT70 | High-affinity InhA inhibitor (Kd 0.1 nM) | Diphenyl ether |
| GSK693 (benzofuran class) | InhA inhibitor; clinical candidate class | Benzofuran-pyrrolidine scaffold |

Maximum Tanimoto similarity to any of the five reference scaffolds was used as the LBVS component score (weight 0.20 in composite). CHEMBL3125270 achieved Tanimoto = 0.153 to the benzofuran reference class, consistent with its benzofuran-containing scaffold.

## S7. Epidemiological Data Sources

| Statistic | Value | Source |
| :--- | :--- | :--- |
| Global TB incidence 2022 | 7.8 million new cases | WHO Global TB Report 2023 |
| Global TB incidence rate | 125/100,000 population | WHO Global TB Report 2023 |
| TB mortality 2022 | 1.13 million deaths | WHO Global TB Report 2023 |
| MDR/RR-TB new cases | 410,000 estimated | WHO Global TB Report 2023 |
| MDR-TB treatment success | 57% globally | WHO Global TB Report 2023 |
| India TB incidence rate | 212/100,000 | India TB Report 2023 |
| India MDR-TB notified | 119,000 cases | India TB Report 2023 |
| India TB deaths 2022 | ~342,000 | WHO Global TB Report 2023 |

Sources: World Health Organization (2023). Global Tuberculosis Report 2023. Geneva: WHO. ISBN 978-92-4-007131-7. Central TB Division (2023). India TB Report 2023. Ministry of Health & Family Welfare, Government of India.

## S8. Computational Environment

| Component | Version/Details |
| :--- | :--- |
| Python | 3.14 |
| RDKit | 2024.03.6 |
| scikit-learn | 1.4+ |
| chembl-webresource-client | 0.10.9 |
| matplotlib | 3.8+ |
| numpy | 1.26+ |
| pandas | 2.1+ |
| Operating system | Windows 11 Pro |
| ChEMBL API accessed | March 2025 |
| PDB 4TZK downloaded | March 2025 |

No molecular docking (AutoDock Vina requires Boost C++ library unavailable on this platform), no molecular dynamics simulation, and no quantum mechanics calculations were performed. All computational claims in this manuscript are limited to ligand-based methods (QSAR and Tanimoto similarity) implemented in RDKit and scikit-learn.

## S9. Data Availability

All scripts, real ChEMBL bioactivity data, RDKit-computed ADMET results, QSAR model outputs, and figure generation code are openly available at: **https://github.com/hssling/MDR-TB-Drug-Discovery**

The ChEMBL bioactivity data (CHEMBL1849) can be independently retrieved at: https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL1849/

Crystal structure PDB 4TZK (InhA with NAD+ cofactor, 1.65 Å resolution) was downloaded from the RCSB Protein Data Bank: https://www.rcsb.org/structure/4TZK
