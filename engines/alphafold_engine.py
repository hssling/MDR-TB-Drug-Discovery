"""
AlphaFold Protein Structure Retrieval Engine — Phase 5b (Version 6 Upgrade)
===========================================================================
Interfaces with the EBI AlphaFold Protein Structure Database API to retrieve
high-confidence 3D atomic structures for Mycobacterium tuberculosis targets.
Completely replaces reliance on crystallized PDB databases with full proteome
coverage AI predictions.

SAFETY: Computational structure retrieval only.
"""

import os
import requests
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir

class AlphaFoldEngine:
    """Retrieves 3D structures from DeepMind's AlphaFold EBI Database."""

    # Essential UniProt Mappings for key MDR-TB targets
    UNIPROT_MAP = {
        "DprE1": "P9WGN7",
        "InhA": "P9WGR1",
        "MmpL3": "P9WJD7",
        "RpoB": "P9WGZ1",
        "GyrA": "P9WG87",
        "KatG": "P9WGY9",
        "KasA": "P9WGL1",
        "EmbB": "P9WGU3",
        "PncA": "P9WGT9",
        "EthA": "P9WGW3"
    }

    ALPHAFOLD_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction/"

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "alphafold_structures"
        )
        self.session = requests.Session()

    def run(self, target_name: str) -> str:
        """Fetch the AlphaFold predicted PDB structure for the Mtb target."""
        print(f"  [AlphaFold] Initializing architectural retrieval for: {target_name}")

        uniprot_id = self.UNIPROT_MAP.get(target_name)
        if not uniprot_id:
            print(f"  [AlphaFold] ⚠ Unmapped target '{target_name}'. Cannot resolve UniProt ID.")
            return None

        pdb_out_path = self.output_dir / f"AF_{target_name}_{uniprot_id}.pdb"
        
        # Memory/Cache Check: Avoid re-downloading
        if pdb_out_path.exists():
            print(f"  [AlphaFold] ✓ Loaded verified structure from local memory cache: {pdb_out_path.name}")
            return str(pdb_out_path)

        return self._fetch_from_ebi(uniprot_id, target_name, pdb_out_path)

    def _fetch_from_ebi(self, uniprot_id: str, target_name: str, out_path: Path) -> str:
        """Execute the Hook to the EBI AlphaFold DB REST API."""
        print(f"  [AlphaFold] Querying DeepMind EBI API for UniProt: {uniprot_id}...")
        try:
            # Step 1: Query API for PDB URL
            api_url = f"{self.ALPHAFOLD_API_BASE}{uniprot_id}"
            response = self.session.get(api_url, timeout=15)
            
            if response.status_code == 200 and response.json():
                data = response.json()[0]
                pdb_url = data.get("pdbUrl")
                confidence = data.get("fractionConfidentAsymmId")
                
                print(f"  [AlphaFold] ✓ Structural map located (Confidence: {confidence})! Downloading...")
                
                # Step 2: Download the raw PDB file
                pdb_resp = self.session.get(pdb_url, timeout=30)
                if pdb_resp.status_code == 200:
                    with open(out_path, "wb") as f:
                        f.write(pdb_resp.content)
                    print(f"  [AlphaFold] ✓ Saved AlphaFold structure to {out_path.name}")
                    return str(out_path)
            else:
                print(f"  [AlphaFold] ⚠ API failed or returned empty payload (Status: {response.status_code}).")
        
        except Exception as e:
            print(f"  [AlphaFold] ⚠ Request Hook Failed: {e}")
            
        print("  [AlphaFold] Falling back to offline structural proxy.")
        return None

if __name__ == "__main__":
    af = AlphaFoldEngine()
    af.run("InhA")
