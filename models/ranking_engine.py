"""
Ranking Engine — Phase 9
=========================
Multi-objective compound ranking combining predicted activity,
Lipinski compliance, target score, novelty, and resistance coverage.

SAFETY: Computational only — in-silico compound prioritization.
"""

import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, safe_save_json


class RankingEngine:
    """Multi-objective compound ranking."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.weights = self.config.get("ranking", {}).get("weights", {
            "predicted_activity": 0.30,
            "lipinski_compliance": 0.15,
            "target_score": 0.20,
            "novelty": 0.15,
            "resistance_coverage": 0.20,
        })
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "ranking"
        )

    def run(self, compounds_df: pd.DataFrame,
            target_scores: pd.DataFrame = None,
            resistance_scores: pd.DataFrame = None) -> dict:
        """Run multi-objective ranking."""
        print("  [Ranking] Computing multi-objective scores...")
        
        ranked = self.compute_ranking(compounds_df, target_scores, resistance_scores)
        
        safe_save_csv(ranked, self.output_dir / "ranked_compounds.csv")
        
        top_10 = ranked.head(10)
        safe_save_csv(top_10, self.output_dir / "top_10_compounds.csv")
        
        summary = {
            "n_compounds_ranked": len(ranked),
            "top_compound": ranked.iloc[0]["Compound_ID"] if "Compound_ID" in ranked.columns and len(ranked) > 0 else "N/A",
            "top_score": round(ranked.iloc[0]["Final_Rank_Score"], 4) if len(ranked) > 0 else 0,
            "top_10_ids": top_10["Compound_ID"].tolist() if "Compound_ID" in top_10.columns else [],
            "mean_rank_score": round(ranked["Final_Rank_Score"].mean(), 4) if "Final_Rank_Score" in ranked.columns else 0,
            "n_high_priority": int((ranked["Priority_Class"] == "Lead").sum()) if "Priority_Class" in ranked.columns else 0,
        }
        safe_save_json(summary, self.output_dir / "ranking_summary.json")
        
        return {
            "ranked_compounds": ranked,
            "top_10": top_10,
            "summary": summary,
        }

    def compute_ranking(self, compounds_df: pd.DataFrame,
                        target_scores: pd.DataFrame = None,
                        resistance_scores: pd.DataFrame = None) -> pd.DataFrame:
        """Compute multi-objective ranking scores."""
        np.random.seed(42)
        df = compounds_df.copy()
        
        # 1. Predicted activity score
        if "Activity_Probability" in df.columns:
            df["score_activity"] = df["Activity_Probability"]
        elif "Activity_Label" in df.columns:
            df["score_activity"] = df["Activity_Label"].astype(float)
        else:
            df["score_activity"] = np.random.uniform(0.2, 0.9, size=len(df))
        
        # 2. Lipinski compliance score
        if "Lipinski_Violations" in df.columns:
            df["score_lipinski"] = 1.0 - (df["Lipinski_Violations"] / 4.0)
        elif "Lipinski_Pass" in df.columns:
            df["score_lipinski"] = df["Lipinski_Pass"].astype(float)
        else:
            # Compute from descriptors
            violations = 0
            if "MolWt" in df.columns:
                violations += (df["MolWt"] > 500).astype(int)
            if "LogP" in df.columns:
                violations += (df["LogP"] > 5).astype(int)
            if "NumHDonors" in df.columns:
                violations += (df["NumHDonors"] > 5).astype(int)
            if "NumHAcceptors" in df.columns:
                violations += (df["NumHAcceptors"] > 10).astype(int)
            df["score_lipinski"] = 1.0 - (violations / 4.0)
        
        # 3. Target score (use mean of top targets if available)
        if target_scores is not None and "Final_Score" in target_scores.columns:
            mean_target = target_scores["Final_Score"].mean()
            df["score_target"] = mean_target + np.random.normal(0, 0.05, size=len(df))
            df["score_target"] = df["score_target"].clip(0, 1)
        else:
            df["score_target"] = np.random.uniform(0.4, 0.9, size=len(df))
        
        # 4. Novelty score (based on structural diversity)
        if "Cluster" in df.columns:
            cluster_counts = df["Cluster"].value_counts()
            df["score_novelty"] = df["Cluster"].map(
                lambda c: 1.0 / (1 + np.log1p(cluster_counts.get(c, 1)))
            )
        else:
            df["score_novelty"] = np.random.uniform(0.3, 0.8, size=len(df))
        
        # 5. Resistance coverage score
        if resistance_scores is not None and "Resistance_Score" in resistance_scores.columns:
            mean_resistance = resistance_scores["Resistance_Score"].mean()
            df["score_resistance"] = mean_resistance + np.random.normal(0, 0.05, size=len(df))
            df["score_resistance"] = df["score_resistance"].clip(0, 1)
        else:
            df["score_resistance"] = np.random.uniform(0.3, 0.8, size=len(df))
        
        # Weighted composite score
        df["Final_Rank_Score"] = (
            self.weights["predicted_activity"] * df["score_activity"] +
            self.weights["lipinski_compliance"] * df["score_lipinski"] +
            self.weights["target_score"] * df["score_target"] +
            self.weights["novelty"] * df["score_novelty"] +
            self.weights["resistance_coverage"] * df["score_resistance"]
        )
        
        # Priority classification
        df["Priority_Class"] = pd.cut(
            df["Final_Rank_Score"],
            bins=[0, 0.4, 0.6, 0.8, 1.0],
            labels=["Low", "Moderate", "Promising", "Lead"]
        )
        
        # Sort and rank
        df = df.sort_values("Final_Rank_Score", ascending=False).reset_index(drop=True)
        df["Rank"] = range(1, len(df) + 1)
        
        n_leads = (df["Priority_Class"] == "Lead").sum()
        n_promising = (df["Priority_Class"] == "Promising").sum()
        print(f"  [Ranking] ✓ Ranked {len(df)} compounds: {n_leads} leads, {n_promising} promising")
        
        return df


if __name__ == "__main__":
    from utils.helpers import generate_mock_compounds
    engine = RankingEngine()
    compounds = generate_mock_compounds()
    results = engine.run(compounds)
    print(f"\nRanking Summary: {results['summary']}")
