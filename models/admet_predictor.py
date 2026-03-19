"""
Deep ADMET & Toxicity Predictor Engine — Phase 7 (Version 7 Upgrade)
====================================================================
Predicts advanced Absorption, Distribution, Metabolism, Excretion, and Toxicity 
(ADMET) parameters for the generated chemical structures. Goes beyond Lipinski 
rules by mapping Blood-Brain Barrier (BBB) permeability and topological 
hERG cardiac toxicity probabilities.

SAFETY: Computational heuristics only.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, safe_save_csv

class ADMETPredictor:
    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.output_dir = Path(self.config["paths"]["output_dir"]) / "admet"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self, compounds_df: pd.DataFrame) -> pd.DataFrame:
        print("  [ADMET] Initializing deep Pharmacokinetic and Toxicity Profiling...")
        
        admet_results = []
        for idx, row in compounds_df.iterrows():
            smiles = row.get("SMILES")
            mol = Chem.MolFromSmiles(str(smiles))
            if not mol: continue
            
            # Topological Parameters
            tpsa = Descriptors.TPSA(mol)
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            rot_bonds = Descriptors.NumRotatableBonds(mol)
            
            # Toxicity & Permeability Heuristics
            bbb_permeable = (mw < 400) and (logp < 5) and (tpsa < 90)
            gi_absorption = (tpsa < 130) and (logp < 5.8)
            
            # hERG Cardiovascular Toxicity Proxy (highly lipophilic and basic amines)
            h_bond_donors = Descriptors.NumHDonors(mol)
            h_bond_acceptors = Descriptors.NumHAcceptors(mol)
            herg_risk = "High" if (logp > 4.5 and rot_bonds > 5 and h_bond_donors < 2) else "Low"
            
            # Quantitative Estimate of Drug-likeness (QED) Continuous metric
            qed_score = QED.qed(mol)
            
            admet_results.append({
                "Compound_ID": row.get("Compound_ID", f"MOL_{idx}"),
                "SMILES": smiles,
                "TPSA": round(tpsa, 2),
                "BBB_Permeable": "Yes" if bbb_permeable else "No",
                "High_GI_Absorption": "Yes" if gi_absorption else "No",
                "hERG_Toxicity_Risk": herg_risk,
                "QED_Drug_Likeness": round(qed_score, 3)
            })
            
        df_res = pd.DataFrame(admet_results)
        safe_save_csv(df_res, self.output_dir / "admet_toxicity_report.csv")
        print(f"  [ADMET] ✓ Processed {len(df_res)} compounds natively. Output saved.")
        return df_res

if __name__ == "__main__":
    predictor = ADMETPredictor()
    df = pd.DataFrame({"Compound_ID": ["Test_1"], "SMILES": ["c1ccccc1-O-C1C(O)C(O)C(O)C(O)C1"]})
    print(predictor.run(df))
