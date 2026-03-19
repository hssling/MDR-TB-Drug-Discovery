"""
DrugBank Data Loader
=====================
Loads drug information using DrugBank-style schema.
Falls back to mock entries for offline operation.

SAFETY: Computational only — uses public drug metadata.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import (
    load_config, ensure_dir, safe_save_csv, safe_save_json,
    generate_mock_drugbank_entries
)


class DrugBankLoader:
    """Load and parse DrugBank-style drug entries."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.db_cfg = self.config.get("data_connectors", {}).get("drugbank", {})
        self.schema_file = self.db_cfg.get("schema_file", "data/mock/drugbank_schema.json")
        self.fallback = self.db_cfg.get("fallback_to_mock", True)
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["data_dir"]) / "drugbank"
        )

    def load(self) -> pd.DataFrame:
        """Load DrugBank entries from schema file or generate mock data."""
        print("  [DrugBank] Loading drug entries...")

        # Try loading from schema file
        if os.path.exists(self.schema_file):
            try:
                with open(self.schema_file, "r") as f:
                    entries = json.load(f)
                df = pd.DataFrame(entries)
                print(f"  [DrugBank] ✓ Loaded {len(df)} drug entries")
                safe_save_csv(df, self.output_dir / "drugbank_entries.csv")
                return df
            except Exception as e:
                print(f"  [DrugBank] ⚠ Error reading schema: {e}")

        # Fallback to mock
        if self.fallback:
            print("  [DrugBank] → Generating mock DrugBank entries")
            entries = generate_mock_drugbank_entries(n=30)
            df = pd.DataFrame(entries)
            safe_save_csv(df, self.output_dir / "drugbank_entries.csv")
            safe_save_json(entries, self.output_dir / "drugbank_schema.json")
            return df
        else:
            raise RuntimeError("DrugBank data unavailable and fallback disabled.")

    def get_tb_drugs(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter for approved TB drugs only."""
        if "status" in df.columns:
            return df[df["status"] == "approved"].copy()
        return df

    def get_targets(self, df: pd.DataFrame) -> list:
        """Extract unique target genes from drug entries."""
        if "target_gene" in df.columns:
            return df["target_gene"].unique().tolist()
        return []


if __name__ == "__main__":
    loader = DrugBankLoader()
    df = loader.load()
    print(f"\nDrugBank entries: {df.shape}")
    print(f"TB drugs: {len(loader.get_tb_drugs(df))}")
    print(f"Targets: {loader.get_targets(df)}")
