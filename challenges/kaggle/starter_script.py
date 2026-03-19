"""
Kaggle Starter Code: MDR-TB Binding Affinity Prediction
=======================================================
This script demonstrates how to parse the pipeline dataset directly inside a Kaggle Notebook,
featurize SMILES using RDKit, and train a baseline Random Forest Regressor to predict
AutoDock Vina thermodynamic scoring (kcal/mol) without running structural simulations.
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# Step 1: Load the docking results derived against AlphaFold InhA
df = pd.read_csv('../../outputs/docking/docking_results_InhA.csv').dropna()
print(f"Dataset Size: {len(df)} molecules")

# Step 2: Extract basic 1D topological RDKit Features (Baseline)
def extract_rdkit_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return [np.nan] * 4
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol)
    ]

# Apply extracting logic
features = df['SMILES'].apply(extract_rdkit_features)
X = pd.DataFrame(features.tolist(), columns=["MolWt", "LogP", "HBD", "HBA"])
y = df['Binding_Affinity_kcal_mol']

# Dropping generation parsing errors 
X = X.dropna()
y = y[X.index]

# Step 3: Train-Test Split (Target predicting simulated docking binding energy)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Step 4: Baseline Regressor
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# Step 5: Evaluate 
predictions = model.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, predictions))

print(f"\nBaseline RMSE (Random Forest on RDKit 1D descr.): {rmse:.4f} kcal/mol")
print("\n--- NEXT STEPS FOR KAGGLERS ---")
print("1. Replace RDKit features with DeepChem/ChemBERTa embeddings.")
print("2. Construct a Graph Neural Network (GNN) using PyTorch Geometric.")
print("3. Try to hit an RMSE of < 0.5 kcal/mol!")
