"""
WHO TB Dataset Parser
======================
Parses WHO tuberculosis surveillance data.
Falls back to synthetic mock data for offline operation.

SAFETY: Computational only — uses public health statistics.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, generate_mock_who_data


class WHOParser:
    """Parse WHO TB incidence/mortality/MDR data."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.who_cfg = self.config.get("data_connectors", {}).get("who", {})
        self.data_file = self.who_cfg.get("data_file", "data/mock/who_tb_data.csv")
        self.fallback = self.who_cfg.get("fallback_to_mock", True)
        self.epi_cfg = self.config.get("epidemiology", {})
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["data_dir"]) / "who"
        )

    def load(self) -> pd.DataFrame:
        """Load WHO TB data from file or generate mock data.
        
        Uses Indian state-level mock data by default for consistency
        with the epidemiology engine's expected schema.
        """
        print("  [WHO] Loading TB surveillance data...")

        # Try loading from existing processed file with expected schema
        if os.path.exists(self.data_file):
            try:
                df = pd.read_csv(self.data_file)
                # Verify it has our expected columns
                if "Region" in df.columns and "Year" in df.columns:
                    print(f"  [WHO] ✓ Loaded {len(df)} records from {self.data_file}")
                    return df
                else:
                    print(f"  [WHO] ⚠ File exists but has incompatible schema, regenerating")
            except Exception as e:
                print(f"  [WHO] ⚠ Error reading {self.data_file}: {e}")

        # Check for previously saved compatible data
        cached = self.output_dir / "who_tb_data.csv"
        if cached.exists():
            try:
                df = pd.read_csv(cached)
                if "Region" in df.columns:
                    print(f"  [WHO] ✓ Loaded {len(df)} cached records")
                    return df
            except Exception:
                pass

        # Generate Indian state-level mock data (compatible with epi engine)
        regions = self.epi_cfg.get("regions", None)
        years = self.epi_cfg.get("years", None)
        print(f"  [WHO] → Generating Indian state-level TB data for {len(regions or [])} regions")
        df = generate_mock_who_data(regions=regions, years=years)
        safe_save_csv(df, self.output_dir / "who_tb_data.csv")
        return df

    def _fetch_who_api(self) -> pd.DataFrame:
        """Attempt to fetch data from WHO Global TB Programme API."""
        try:
            import requests
            url = "https://extranet.who.int/tme/generateCSV.asp?ds=notifications"
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                from io import StringIO
                df = pd.read_csv(StringIO(resp.text))
                print(f"  [WHO] ✓ Fetched {len(df)} records from WHO API")
                return df
        except Exception:
            pass
        return None

    def get_summary(self, df: pd.DataFrame) -> dict:
        """Generate summary statistics from WHO data."""
        summary = {
            "total_records": len(df),
            "regions": df["Region"].nunique() if "Region" in df.columns else 0,
            "year_range": f"{df['Year'].min()}-{df['Year'].max()}" if "Year" in df.columns else "N/A",
            "mean_incidence": df["TB_Incidence_per_100k"].mean() if "TB_Incidence_per_100k" in df.columns else 0,
            "mean_mdr_pct": df["MDR_TB_Percentage"].mean() if "MDR_TB_Percentage" in df.columns else 0,
        }
        return summary


if __name__ == "__main__":
    parser = WHOParser()
    df = parser.load()
    print(f"\nWHO data shape: {df.shape}")
    print(parser.get_summary(df))
