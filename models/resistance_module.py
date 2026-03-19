"""
Resistance Module — Phase 8
=============================
MDR gene mapping, mutation cataloging, and resistance scoring.

SAFETY: Computational only — bioinformatic analysis of resistance markers.
"""

import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, safe_save_json


class ResistanceModule:
    """Map MDR-TB resistance genes and compute resistance scores."""

    # Comprehensive MDR-TB mutation catalog
    MUTATION_CATALOG = {
        "rpoB": {
            "mutations": ["S450L", "H445Y", "D435V", "S450W", "L452P", "H445D"],
            "drug": "Rifampicin",
            "frequency": [0.45, 0.12, 0.08, 0.05, 0.04, 0.03],
            "confidence": "high",
        },
        "katG": {
            "mutations": ["S315T", "S315N", "S315I"],
            "drug": "Isoniazid",
            "frequency": [0.65, 0.05, 0.02],
            "confidence": "high",
        },
        "inhA": {
            "mutations": ["C-15T", "T-8A", "T-8C", "S94A"],
            "drug": "Isoniazid (low-level)",
            "frequency": [0.25, 0.05, 0.03, 0.04],
            "confidence": "high",
        },
        "embB": {
            "mutations": ["M306V", "M306I", "G406A", "Q497R"],
            "drug": "Ethambutol",
            "frequency": [0.35, 0.15, 0.08, 0.06],
            "confidence": "moderate",
        },
        "pncA": {
            "mutations": ["D49N", "H51D", "L85R", "V139A", "T135P"],
            "drug": "Pyrazinamide",
            "frequency": [0.08, 0.06, 0.05, 0.04, 0.03],
            "confidence": "moderate",
        },
        "gyrA": {
            "mutations": ["D94G", "A90V", "D94A", "D94N", "D94Y"],
            "drug": "Fluoroquinolones",
            "frequency": [0.30, 0.20, 0.10, 0.08, 0.05],
            "confidence": "high",
        },
        "rrs": {
            "mutations": ["A1401G", "C1402T", "G1484T"],
            "drug": "Aminoglycosides",
            "frequency": [0.50, 0.10, 0.08],
            "confidence": "high",
        },
        "ethA": {
            "mutations": ["A381P", "Q254Stop", "Y84D"],
            "drug": "Ethionamide",
            "frequency": [0.12, 0.08, 0.05],
            "confidence": "moderate",
        },
        "Rv0678": {
            "mutations": ["V1F", "M1R", "R109Stop"],
            "drug": "Bedaquiline/Clofazimine",
            "frequency": [0.10, 0.08, 0.05],
            "confidence": "moderate",
        },
        "atpE": {
            "mutations": ["E61D", "A63P", "I66M", "D28V"],
            "drug": "Bedaquiline",
            "frequency": [0.15, 0.10, 0.08, 0.05],
            "confidence": "high",
        },
    }

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.res_cfg = self.config.get("resistance", {})
        self.mdr_genes_config = self.res_cfg.get("mdr_genes", [])
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "resistance"
        )

    def run(self, de_results: pd.DataFrame = None) -> dict:
        """Run full resistance analysis."""
        print("  [Resistance] Mapping MDR genes...")
        gene_map = self.map_mdr_genes()
        
        print("  [Resistance] Building mutation catalog...")
        mutations = self.build_mutation_catalog()
        
        print("  [Resistance] Computing resistance scores...")
        scores = self.compute_resistance_scores(de_results)
        
        # Save results
        safe_save_csv(gene_map, self.output_dir / "mdr_gene_map.csv")
        safe_save_csv(mutations, self.output_dir / "mutation_catalog.csv")
        safe_save_csv(scores, self.output_dir / "resistance_scores.csv")
        
        summary = {
            "n_resistance_genes": len(gene_map),
            "n_mutations_cataloged": len(mutations),
            "drugs_covered": gene_map["Drug"].nunique() if "Drug" in gene_map.columns else 0,
            "high_confidence_genes": int((gene_map["Confidence"] == "high").sum()) if "Confidence" in gene_map.columns else 0,
            "average_resistance_score": round(scores["Resistance_Score"].mean(), 4) if "Resistance_Score" in scores.columns else 0,
        }
        safe_save_json(summary, self.output_dir / "resistance_summary.json")
        
        return {
            "gene_map": gene_map,
            "mutations": mutations,
            "scores": scores,
            "summary": summary,
        }

    def map_mdr_genes(self) -> pd.DataFrame:
        """Map MDR resistance genes with their associated drugs and functions."""
        rows = []
        
        for gene_name, info in self.MUTATION_CATALOG.items():
            # Check for config weight
            config_weight = 0.5
            for cfg_gene in self.mdr_genes_config:
                if cfg_gene.get("gene", "").lower() == gene_name.lower():
                    config_weight = cfg_gene.get("weight", 0.5)
                    break
            
            rows.append({
                "Gene": gene_name,
                "Drug": info["drug"],
                "N_Known_Mutations": len(info["mutations"]),
                "Top_Mutation": info["mutations"][0] if info["mutations"] else "N/A",
                "Top_Mutation_Frequency": info["frequency"][0] if info["frequency"] else 0,
                "Confidence": info["confidence"],
                "Config_Weight": config_weight,
            })
        
        df = pd.DataFrame(rows)
        print(f"  [Resistance] ✓ Mapped {len(df)} resistance genes")
        return df

    def build_mutation_catalog(self) -> pd.DataFrame:
        """Build comprehensive mutation catalog."""
        rows = []
        for gene_name, info in self.MUTATION_CATALOG.items():
            for mut, freq in zip(info["mutations"], info["frequency"]):
                rows.append({
                    "Gene": gene_name,
                    "Mutation": mut,
                    "Drug": info["drug"],
                    "Frequency": freq,
                    "Confidence": info["confidence"],
                    "Clinical_Significance": "Resistance-associated" if freq > 0.1 else "Possible resistance",
                })
        
        df = pd.DataFrame(rows)
        df = df.sort_values(["Gene", "Frequency"], ascending=[True, False]).reset_index(drop=True)
        print(f"  [Resistance] ✓ Cataloged {len(df)} mutations across {df['Gene'].nunique()} genes")
        return df

    def compute_resistance_scores(self, de_results: pd.DataFrame = None) -> pd.DataFrame:
        """Compute resistance score for each gene combining multiple factors."""
        np.random.seed(42)
        rows = []
        
        for gene_name, info in self.MUTATION_CATALOG.items():
            n_mutations = len(info["mutations"])
            max_freq = max(info["frequency"]) if info["frequency"] else 0
            confidence_score = 1.0 if info["confidence"] == "high" else 0.7
            
            # Expression association (from DE results if available)
            expr_score = 0.5
            if de_results is not None and "Gene" in de_results.columns:
                matched = de_results[de_results["Gene"].astype(str).str.contains(gene_name[:4], case=False, na=False)]
                if len(matched) > 0:
                    fc = abs(matched.iloc[0].get("log2FC", 0))
                    expr_score = min(fc / 3.0, 1.0)
            
            # Composite resistance score
            res_score = (
                0.30 * max_freq +
                0.25 * (n_mutations / 10) +
                0.25 * confidence_score +
                0.20 * expr_score
            )
            
            rows.append({
                "Gene": gene_name,
                "Drug": info["drug"],
                "N_Mutations": n_mutations,
                "Max_Mutation_Frequency": round(max_freq, 3),
                "Confidence_Score": round(confidence_score, 3),
                "Expression_Score": round(expr_score, 3),
                "Resistance_Score": round(res_score, 4),
                "Priority": "Critical" if res_score > 0.5 else ("High" if res_score > 0.3 else "Moderate"),
            })
        
        df = pd.DataFrame(rows)
        df = df.sort_values("Resistance_Score", ascending=False).reset_index(drop=True)
        print(f"  [Resistance] ✓ Scored {len(df)} resistance genes")
        return df


if __name__ == "__main__":
    module = ResistanceModule()
    results = module.run()
    print(f"\nResistance Summary: {results['summary']}")
