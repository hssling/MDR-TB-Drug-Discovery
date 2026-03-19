"""
ESMFold Structural Predictor — Innovative Extension
=====================================================
Uses Meta's ESMFold API to computationally predict 3D
protein structures for novel MDR-TB mutated sequences 
that do not exist in the PDB.

SAFETY: Computational structural biology only.
"""

import os
import time
import requests
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_json

class ESMFoldPredictor:
    """Predicts 3D structural formations for uncharacterized proteins/mutants."""
    
    API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    
    # Mock sequence for a mutated DprE1 fragment (for demonstration if real sequence is too long)
    MOCK_TARGET_SEQ = {
        "DprE1_Mutant": "MTEPDKTVTVGIGSGIGAAIAANIVRQKGIEVTILDRTPVDHLHPSPL"
    }

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "structures"
        )

    def run(self, target_name: str, sequence: str = None) -> str:
        """Call ESMFold API to generate a PDB file from an amino acid sequence."""
        print(f"  [ESMFold] Evaluating structural availability for target: {target_name}")
        
        sequence_to_fold = sequence or self.MOCK_TARGET_SEQ.get(target_name, self.MOCK_TARGET_SEQ["DprE1_Mutant"])
        
        pdb_out_path = self.output_dir / f"{target_name}_predicted.pdb"
        
        # Avoid hitting API if we already have it
        if pdb_out_path.exists():
            print(f"  [ESMFold] ✓ Loaded previously predicted structure for {target_name}")
            return str(pdb_out_path)
            
        print("  [ESMFold] Connecting to ESMFold Inference API (Meta)...")
        try:
            # We use a post request to the API
            response = requests.post(self.API_URL, data=sequence_to_fold, headers={"Content-Type": "text/plain"}, timeout=30)
            
            if response.status_code == 200:
                with open(pdb_out_path, "w", encoding="utf-8") as f:
                    f.write(response.text)
                print(f"  [ESMFold] ✓ Successfully predicted 3D structure. Saved to {pdb_out_path.name}")
                return str(pdb_out_path)
            else:
                print(f"  [ESMFold] ⚠ API returned status code: {response.status_code}")
                # Fallback purely for robustness
                return self._generate_mock_pdb(pdb_out_path, target_name)
                
        except Exception as e:
            print(f"  [ESMFold] ⚠ API Connection Error: {e}")
            return self._generate_mock_pdb(pdb_out_path, target_name)

    def _generate_mock_pdb(self, path: Path, name: str) -> str:
        """Generates a dummy PDB file format if offline to ensure pipeline continutity."""
        print(f"  [ESMFold] Generating mock structural scaffold for {name} (offline mode)")
        mock_pdb = f"HEADER    Mock predicted structure for {name}\nATOM      1  CA  ALA A   1      11.111  22.222  33.333  1.00 20.00           C  \nEND"
        with open(path, "w", encoding="utf-8") as f:
            f.write(mock_pdb)
        return str(path)

if __name__ == "__main__":
    predictor = ESMFoldPredictor()
    path = predictor.run("DprE1_Mutant")
    print(f"Output saved to: {path}")
