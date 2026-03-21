# MDR-TB InhA Inhibitor Ranking Challenge

**Kaggle Dataset**: [jkhospital/mdrtb-drug-discovery-ai-pipeline](https://www.kaggle.com/datasets/jkhospital/mdrtb-drug-discovery-ai-pipeline)

---

## Overview

Multidrug-resistant tuberculosis (MDR-TB) killed an estimated 1.13 million people in 2023 (WHO Global TB Report). The enoyl-ACP reductase enzyme **InhA** (UniProt: P9WGR1, ChEMBL target: CHEMBL1849) is a clinically validated first-line drug target — inhibited by isoniazid after activation by KatG.

This challenge provides a **real, experimentally measured dataset** of 421 InhA inhibitors from ChEMBL with their IC50 values, physicochemical properties, and composite LBVS ranking scores. Your goal is to **predict which compounds will rank highest** and understand what molecular features drive potency.

---

## The Dataset

### Primary File: `outputs/ranking/ranked_compounds.csv`

421 compounds with measured InhA IC50 values from ChEMBL1849, plus computed properties.

| Column | Type | Description |
|--------|------|-------------|
| `molecule_chembl_id` | string | ChEMBL compound ID |
| `canonical_smiles` | string | RDKit-canonical SMILES |
| `IC50_nM` | float | **Measured** IC50 in nM (from published ChEMBL assays) |
| `pIC50` | float | −log₁₀(IC50 in M) — higher = more potent |
| `MW` | float | Molecular weight (Da) |
| `LogP` | float | Calculated lipophilicity |
| `TPSA` | float | Topological polar surface area (Ų) |
| `HBD` | int | Hydrogen bond donors |
| `HBA` | int | Hydrogen bond acceptors |
| `RotBonds` | int | Rotatable bonds |
| `Rings` | int | Ring count |
| `QED` | float | Drug-likeness score [0–1] (Bickerton 2012) |
| `Lipinski_Pass` | bool | Lipinski Rule-of-Five compliance |
| `hERG_Risk` | str | Cardiotoxicity flag (Low/High) |
| `GI_Absorption` | str | GI absorption proxy |
| `LBVS_Composite_Score` | float | Ligand-based virtual screening score [0–1] |
| `QSAR_Active_Prob` | float | ML-predicted activity probability [0–1] |
| `Composite_Score` | float | Final ranking: 0.40×pIC50 + 0.30×QSAR + 0.20×LBVS + 0.10×QED |
| `Rank` | int | Compound rank (1 = best) |

**Key fact**: All IC50 values are real, experimentally measured values from published ChEMBL assays — not computationally predicted.

### Supporting Files

- `outputs/ranking/top_10_compounds.csv` — Top 10 ranked compounds
- `outputs/ranking/ranking_summary.json` — Pipeline run metadata and statistics
- `outputs/manuscript/manuscript_v8_genuine.md` — Full peer-reviewed manuscript describing the methods

---

## Challenge Tasks

### Task 1 (Beginner) — Reproduce the Composite Score

Use the ranking formula to recompute `Composite_Score` from its components and verify against the provided rankings:

```
Composite_Score = 0.40 × pIC50_norm + 0.30 × QSAR_Active_Prob + 0.20 × LBVS_Composite_Score + 0.10 × QED
```

Where `pIC50_norm` is min-max normalised pIC50 across all 421 compounds.

**Deliverable**: A notebook showing your recomputed scores correlate with `Composite_Score` (r > 0.99).

### Task 2 (Intermediate) — QSAR Binary Classifier

Train a binary classifier to predict `QSAR_Predicted_Active` (IC50 ≤ 1000 nM = active) from ECFP4 fingerprints:

1. Generate Morgan fingerprints (radius=2, 2048 bits) from `canonical_smiles`
2. Split 80/20 stratified on the activity label
3. Train Random Forest, Gradient Boosting, or Logistic Regression
4. Report ROC-AUC and compare to the pipeline's 0.979

**Bonus**: Try XGBoost or a graph neural network (PyTorch Geometric).

### Task 3 (Advanced) — Improve the Ranking

Propose and implement a better ranking formula. Ideas:
- Replace LBVS Tanimoto with 3D shape similarity
- Add a selectivity term against human off-targets
- Use ChemBERTa embeddings instead of ECFP4
- Incorporate synthetic accessibility (SAScore)

**Deliverable**: A notebook with your improved ranking and justification for why it would perform better in a prospective screen.

---

## Starter Code

See `challenges/kaggle/starter_script.py` for a working baseline that:
- Loads the dataset
- Generates ECFP4 fingerprints
- Trains a Random Forest classifier
- Reports ROC-AUC

---

## How to Use This Dataset in a Kaggle Notebook

```python
import pandas as pd

# The dataset is mounted at /kaggle/input/mdrtb-drug-discovery-ai-pipeline/
df = pd.read_csv('/kaggle/input/mdrtb-drug-discovery-ai-pipeline/outputs/ranking/ranked_compounds.csv')
print(f"Dataset: {len(df)} compounds")
print(df[['molecule_chembl_id', 'IC50_nM', 'pIC50', 'Composite_Score', 'Rank']].head(10))
```

---

## Sharing Guidelines

- Publish your notebook to the dataset page so others can learn from your approach
- Clearly state whether IC50 values are real measured values or predicted
- If you use additional external data, cite your sources
- Do not claim predicted scores as experimental IC50 values

---

## Scientific Context

- **Target**: InhA (2-trans-enoyl-ACP reductase), NADH-dependent, primary target of isoniazid
- **Resistance mechanism**: katG S315T mutation (~67% of MDR strains, CRyPTIC n=12,289) prevents pro-drug activation; direct InhA inhibitors bypass this
- **Reference structure**: PDB 4TZK (InhA with triclosan-class inhibitor, resolution 2.0 Å)
- **Top compound**: CHEMBL3125270 — IC50 4 nM, Composite Score 0.825, Lipinski-compliant

---

*Dataset curated by Dr Siddalingaiah H S, Shridevi Institute of Medical Sciences and Research Hospital, Tumkur. Source code: [github.com/hssling/MDR-TB-Drug-Discovery](https://github.com/hssling/MDR-TB-Drug-Discovery)*
