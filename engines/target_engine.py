"""
Target Engine — Phase 5
========================
Multi-factor drug target scoring for MDR-TB.
Combines druggability, essentiality, resistance association,
expression fold change, and conservation scores.

SAFETY: Computational only — in-silico target prioritization.
"""

import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, safe_save_json


class TargetEngine:
    """Multi-factor scoring of potential drug targets."""

    # Known TB drug targets with properties
    KNOWN_TARGETS = {
        "InhA": {"druggability": 0.92, "essentiality": 0.95, "conservation": 0.98, "category": "Cell wall"},
        "RpoB": {"druggability": 0.88, "essentiality": 0.99, "conservation": 0.97, "category": "Transcription"},
        "KatG": {"druggability": 0.75, "essentiality": 0.80, "conservation": 0.90, "category": "Oxidative stress"},
        "EmbB": {"druggability": 0.85, "essentiality": 0.88, "conservation": 0.92, "category": "Cell wall"},
        "PncA": {"druggability": 0.70, "essentiality": 0.72, "conservation": 0.85, "category": "Metabolism"},
        "GyrA": {"druggability": 0.90, "essentiality": 0.97, "conservation": 0.95, "category": "DNA replication"},
        "GyrB": {"druggability": 0.87, "essentiality": 0.96, "conservation": 0.94, "category": "DNA replication"},
        "AtpE": {"druggability": 0.82, "essentiality": 0.91, "conservation": 0.93, "category": "Energy metabolism"},
        "DprE1": {"druggability": 0.95, "essentiality": 0.93, "conservation": 0.96, "category": "Cell wall"},
        "MmpL3": {"druggability": 0.88, "essentiality": 0.94, "conservation": 0.91, "category": "Cell wall"},
        "Rrs": {"druggability": 0.65, "essentiality": 0.98, "conservation": 0.99, "category": "Translation"},
        "ClpP1P2": {"druggability": 0.78, "essentiality": 0.86, "conservation": 0.88, "category": "Protein degradation"},
        "FadD32": {"druggability": 0.80, "essentiality": 0.89, "conservation": 0.90, "category": "Lipid metabolism"},
        "Pks13": {"druggability": 0.83, "essentiality": 0.90, "conservation": 0.92, "category": "Lipid metabolism"},
        "MurB": {"druggability": 0.77, "essentiality": 0.85, "conservation": 0.89, "category": "Cell wall"},
    }

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.weights = self.config.get("target_scoring", {}).get("weights", {
            "druggability": 0.25,
            "essentiality": 0.25,
            "resistance_association": 0.20,
            "expression_fold_change": 0.15,
            "conservation": 0.15,
        })
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "targets"
        )

    def run(self, de_results: pd.DataFrame = None,
            resistance_genes: list = None) -> dict:
        """Run multi-factor target scoring."""
        print("  [Target] Scoring drug targets...")
        
        scored_targets = self.score_targets(de_results, resistance_genes)
        
        safe_save_csv(scored_targets, self.output_dir / "scored_targets.csv")
        
        summary = {
            "n_targets_scored": len(scored_targets),
            "top_target": scored_targets.iloc[0]["Target"] if len(scored_targets) > 0 else "N/A",
            "top_score": round(scored_targets.iloc[0]["Final_Score"], 4) if len(scored_targets) > 0 else 0,
            "high_priority_targets": scored_targets[scored_targets["Priority"] == "High"]["Target"].tolist(),
        }
        safe_save_json(summary, self.output_dir / "target_summary.json")
        
        return {
            "scored_targets": scored_targets,
            "summary": summary,
        }

    def score_targets(self, de_results: pd.DataFrame = None,
                      resistance_genes: list = None) -> pd.DataFrame:
        """Score each target using multiple factors."""
        np.random.seed(42)
        
        # Build resistance gene set
        res_genes = set()
        if resistance_genes:
            for rg in resistance_genes:
                if isinstance(rg, dict):
                    res_genes.add(rg.get("gene", ""))
                else:
                    res_genes.add(str(rg))
        else:
            res_genes = {"rpoB", "katG", "inhA", "embB", "pncA", "gyrA", "rrs", "ethA"}
        
        results = []
        for target_name, props in self.KNOWN_TARGETS.items():
            druggability = props["druggability"]
            essentiality = props["essentiality"]
            conservation = props["conservation"]
            
            # Resistance association
            resistance_assoc = 0.9 if target_name.lower() in {g.lower() for g in res_genes} else np.random.uniform(0.1, 0.5)
            
            # Expression fold change (from DE results if available)
            expr_fc = 0.5
            if de_results is not None and "Gene" in de_results.columns:
                matched = de_results[de_results["Gene"].astype(str).str.contains(target_name[:4], case=False, na=False)]
                if len(matched) > 0:
                    fc = abs(matched.iloc[0].get("log2FC", 0))
                    expr_fc = min(fc / 4.0, 1.0)  # Normalize to [0, 1]
            
            # Weighted composite score
            final_score = (
                self.weights["druggability"] * druggability +
                self.weights["essentiality"] * essentiality +
                self.weights["resistance_association"] * resistance_assoc +
                self.weights["expression_fold_change"] * expr_fc +
                self.weights["conservation"] * conservation
            )
            
            # Priority classification
            if final_score >= 0.80:
                priority = "High"
            elif final_score >= 0.60:
                priority = "Medium"
            else:
                priority = "Low"
            
            results.append({
                "Target": target_name,
                "Category": props["category"],
                "Druggability": round(druggability, 3),
                "Essentiality": round(essentiality, 3),
                "Conservation": round(conservation, 3),
                "Resistance_Association": round(resistance_assoc, 3),
                "Expression_FC_Score": round(expr_fc, 3),
                "Final_Score": round(final_score, 4),
                "Priority": priority,
            })
        
        df = pd.DataFrame(results)
        df = df.sort_values("Final_Score", ascending=False).reset_index(drop=True)
        df["Rank"] = range(1, len(df) + 1)
        
        n_high = (df["Priority"] == "High").sum()
        print(f"  [Target] ✓ Scored {len(df)} targets → {n_high} high-priority")
        return df


if __name__ == "__main__":
    engine = TargetEngine()
    results = engine.run()
    print(f"\nTarget Summary: {results['summary']}")
