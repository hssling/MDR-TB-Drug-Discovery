"""
PubChem Structure Loader
=========================
Loads compound structures and properties from PubChem.
Falls back to mock compound data for offline operation.

SAFETY: Computational only — uses public chemical databases.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import (
    load_config, ensure_dir, safe_save_csv,
    generate_mock_compounds
)


class PubChemLoader:
    """Load compound data from PubChem REST API or mock data."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.pc_cfg = self.config.get("data_connectors", {}).get("pubchem", {})
        self.compound_ids = self.pc_cfg.get("compound_ids", [2244, 5291])
        self.fallback = self.pc_cfg.get("fallback_to_mock", True)
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["data_dir"]) / "pubchem"
        )

    def load(self, n_compounds: int = 50) -> pd.DataFrame:
        """Load compound data from PubChem or generate mock data."""
        print("  [PubChem] Loading compound structures...")

        # Try real PubChem API
        try:
            df = self._fetch_pubchem()
            if df is not None and not df.empty:
                safe_save_csv(df, self.output_dir / "compounds.csv")
                return df
        except Exception as e:
            print(f"  [PubChem] ⚠ PubChem API unavailable: {e}")

        # Fallback to mock
        if self.fallback:
            print(f"  [PubChem] → Generating {n_compounds} mock compounds")
            df = generate_mock_compounds(n=n_compounds)
            safe_save_csv(df, self.output_dir / "compounds.csv")
            return df
        else:
            raise RuntimeError("PubChem data unavailable and fallback disabled.")

    def _fetch_pubchem(self) -> pd.DataFrame:
        """Fetch compound properties from PubChem REST API."""
        try:
            import requests
            base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
            properties = "MolecularFormula,MolecularWeight,CanonicalSMILES,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount"
            
            cids = ",".join(str(c) for c in self.compound_ids)
            url = f"{base_url}/compound/cid/{cids}/property/{properties}/JSON"
            
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                props = data.get("PropertyTable", {}).get("Properties", [])
                if props:
                    df = pd.DataFrame(props)
                    print(f"  [PubChem] ✓ Fetched {len(df)} compounds from PubChem")
                    return df
        except ImportError:
            print("  [PubChem] requests library not available")
        except Exception as e:
            print(f"  [PubChem] API error: {e}")
        return None

    def get_smiles_list(self, df: pd.DataFrame) -> list:
        """Extract SMILES strings from compound DataFrame."""
        smiles_col = None
        for col in ["SMILES", "CanonicalSMILES", "smiles"]:
            if col in df.columns:
                smiles_col = col
                break
        if smiles_col:
            return df[smiles_col].dropna().tolist()
        return []


if __name__ == "__main__":
    loader = PubChemLoader()
    df = loader.load()
    print(f"\nCompound data shape: {df.shape}")
    print(f"SMILES extracted: {len(loader.get_smiles_list(df))}")
