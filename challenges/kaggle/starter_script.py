"""
Kaggle Starter: MDR-TB InhA Inhibitor Ranking Challenge
========================================================
Predicts InhA inhibitor activity (IC50 <= 1000 nM = active) from ECFP4 fingerprints.

Dataset: https://www.kaggle.com/datasets/jkhospital/mdrtb-drug-discovery-ai-pipeline
Source:  https://github.com/hssling/MDR-TB-Drug-Discovery

Run in a Kaggle notebook:
    The dataset is at /kaggle/input/mdrtb-drug-discovery-ai-pipeline/
    All IC50 values are real, measured values from ChEMBL — not predicted.
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.preprocessing import StandardScaler

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# ── 1. Load the dataset ────────────────────────────────────────────────────────

# In a Kaggle notebook:
DATA_PATH = "/kaggle/input/mdrtb-drug-discovery-ai-pipeline/outputs/ranking/ranked_compounds.csv"

# When running locally from the repo root:
# DATA_PATH = "outputs/ranking/ranked_compounds.csv"

df = pd.read_csv(DATA_PATH)
print(f"Dataset loaded: {len(df)} compounds")
print(f"Active compounds (IC50 ≤ 1000 nM): {(df['IC50_nM'] <= 1000).sum()}")
print(f"IC50 range: {df['IC50_nM'].min():.1f} – {df['IC50_nM'].max():.0f} nM")
print(f"\nTop 5 compounds by Composite Score:")
print(df[['molecule_chembl_id', 'IC50_nM', 'pIC50', 'Composite_Score', 'Rank']].head(5).to_string(index=False))

# ── 2. Generate ECFP4 Morgan fingerprints ─────────────────────────────────────

def smiles_to_ecfp4(smiles: str, n_bits: int = 2048) -> np.ndarray | None:
    """Convert a SMILES string to a 2048-bit ECFP4 Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    return np.array(fp)

print("\nGenerating ECFP4 fingerprints...")
fps = df["canonical_smiles"].apply(smiles_to_ecfp4)
valid_mask = fps.notna()
print(f"Valid SMILES: {valid_mask.sum()} / {len(df)}")

X_fp = np.stack(fps[valid_mask].values)             # ECFP4 features
y = (df.loc[valid_mask, "IC50_nM"] <= 1000).astype(int).values   # binary label

print(f"\nFeature matrix: {X_fp.shape}")
print(f"Active (label=1): {y.sum()}, Inactive (label=0): {(y == 0).sum()}")

# ── 3. Optional: combine fingerprints with RDKit 1D descriptors ───────────────

DESCRIPTOR_COLS = ["MW", "LogP", "TPSA", "HBD", "HBA", "QED", "RotBonds", "Rings"]
X_desc = df.loc[valid_mask, DESCRIPTOR_COLS].fillna(0).values

scaler = StandardScaler()
X_desc_scaled = scaler.fit_transform(X_desc)

X_combined = np.hstack([X_fp, X_desc_scaled])

# ── 4. Train models with 5-fold cross-validation ──────────────────────────────

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

models = {
    "Random Forest (ECFP4)": (RandomForestClassifier(n_estimators=200, random_state=42), X_fp),
    "Random Forest (ECFP4 + descriptors)": (RandomForestClassifier(n_estimators=200, random_state=42), X_combined),
    "Gradient Boosting (ECFP4)": (GradientBoostingClassifier(n_estimators=100, random_state=42), X_fp),
    "Logistic Regression (ECFP4)": (LogisticRegression(max_iter=1000, random_state=42), X_fp),
}

print("\n── 5-Fold Cross-Validation Results ──")
print(f"{'Model':<45} {'ROC-AUC':>10} {'±SD':>8}")
print("-" * 65)

best_model_name = None
best_auc = 0.0
for name, (model, X) in models.items():
    aucs = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
    mean_auc = aucs.mean()
    print(f"{name:<45} {mean_auc:.4f}     ±{aucs.std():.4f}")
    if mean_auc > best_auc:
        best_auc = mean_auc
        best_model_name = name

print(f"\nPipeline baseline ROC-AUC:  0.9790")
print(f"Your best model ({best_model_name}): {best_auc:.4f}")

# ── 5. Full train + evaluation on held-out 20% ────────────────────────────────

from sklearn.model_selection import train_test_split

best_model, best_X = models[best_model_name]
X_train, X_test, y_train, y_test = train_test_split(
    best_X, y, test_size=0.2, stratify=y, random_state=42
)
best_model.fit(X_train, y_train)
y_prob = best_model.predict_proba(X_test)[:, 1]
y_pred = best_model.predict(X_test)

test_auc = roc_auc_score(y_test, y_prob)
print(f"\nHeld-out test ROC-AUC: {test_auc:.4f}")
print("\nClassification report:")
print(classification_report(y_test, y_pred, target_names=["Inactive", "Active"]))

# ── 6. Feature importance (top 20 descriptors) ────────────────────────────────

if hasattr(best_model, "feature_importances_") and best_X is X_combined:
    importances = best_model.feature_importances_
    # Descriptor importances (last N_DESC features)
    desc_importances = importances[-len(DESCRIPTOR_COLS):]
    desc_df = pd.DataFrame({
        "Descriptor": DESCRIPTOR_COLS,
        "Importance": desc_importances
    }).sort_values("Importance", ascending=False)
    print("\nTop descriptor importances:")
    print(desc_df.to_string(index=False))

# ── 7. Next steps ─────────────────────────────────────────────────────────────

print("""
── Suggested Improvements ──────────────────────────────────────────────────────
1. Try ChemBERTa embeddings: from transformers import AutoTokenizer, AutoModel
   → Encode SMILES strings as 600-dim dense vectors (better than ECFP4 for activity)

2. Graph Neural Network: pip install torch-geometric
   → Represent molecules as graphs — captures bond topology better than 2D fingerprints

3. Regression instead of classification: predict pIC50 directly (RMSE metric)

4. Ensemble: combine RF + GBM + LR predictions (voting or stacking)

5. Publish your notebook to the dataset page for the community!
   https://www.kaggle.com/datasets/jkhospital/mdrtb-drug-discovery-ai-pipeline
────────────────────────────────────────────────────────────────────────────────
""")
