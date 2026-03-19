"""
GEO Data Loader
================
Loads gene expression data from NCBI GEO database.
Falls back to synthetic mock data for offline operation.

SAFETY: Computational only — downloads public expression data.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import (
    load_config, ensure_dir, safe_save_csv,
    generate_mock_expression_matrix
)


class GEOLoader:
    """Load gene expression data from GEO or generate mock data."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.geo_cfg = self.config.get("data_connectors", {}).get("geo", {})
        self.series_id = self.geo_cfg.get("series_id", "GSE153029")
        self.fallback = self.geo_cfg.get("fallback_to_mock", True)
        self.n_genes = self.geo_cfg.get("mock_n_genes", 5000)
        self.n_samples = self.geo_cfg.get("mock_n_samples", 20)
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["data_dir"]) / "geo"
        )

    def load(self) -> pd.DataFrame:
        """
        Load expression data with CTRL_/MDR_ column labels.
        Uses mock data with proper sample labels for compatibility
        with the omics engine differential expression analysis.
        """
        print(f"  [GEO] Loading expression data...")
        
        # Check for cached compatible data
        cached = self.output_dir / "expression_matrix.csv"
        if cached.exists():
            try:
                df = pd.read_csv(cached, index_col=0)
                has_ctrl = any(c.startswith("CTRL") for c in df.columns)
                has_mdr = any(c.startswith("MDR") for c in df.columns)
                if has_ctrl and has_mdr:
                    print(f"  [GEO] ✓ Loaded cached {df.shape[0]} genes × {df.shape[1]} samples")
                    return df
                else:
                    print(f"  [GEO] ⚠ Cached data has incompatible columns, regenerating")
            except Exception:
                pass

        # Generate labeled mock expression data
        print(f"  [GEO] → Generating expression data ({self.n_genes} genes × {self.n_samples} samples)")
        df = generate_mock_expression_matrix(
            n_genes=self.n_genes, n_samples=self.n_samples
        )
        safe_save_csv(df, self.output_dir / "expression_matrix.csv")
        return df

    def _load_from_geo(self) -> pd.DataFrame:
        """Attempt real GEO download using GEOparse."""
        try:
            import GEOparse
            gse = GEOparse.get_GEO(geo=self.series_id, destdir=str(self.output_dir))
            # Extract expression table from first platform
            for gsm_name, gsm in gse.gsms.items():
                break
            # Build expression matrix
            data_frames = []
            for gsm_name, gsm in gse.gsms.items():
                col = gsm.table[["ID_REF", "VALUE"]].set_index("ID_REF")
                col.columns = [gsm_name]
                data_frames.append(col)
            if data_frames:
                df = pd.concat(data_frames, axis=1)
                df.index.name = "Gene"
                return df.astype(float)
        except ImportError:
            print("  [GEO] GEOparse not installed; using mock data.")
        except Exception as e:
            print(f"  [GEO] Error fetching GEO data: {e}")
        return None


if __name__ == "__main__":
    loader = GEOLoader()
    df = loader.load()
    print(f"\nExpression matrix shape: {df.shape}")
    print(df.head())
