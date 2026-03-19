"""
Epidemiology Engine — Phase 4
==============================
TB trends analysis, MDR pattern detection, and district-level aggregation.

SAFETY: Computational only — public health statistical modeling.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, safe_save_json


class EpidemiologyEngine:
    """Analyze TB epidemiological trends, MDR patterns, and regional aggregation."""

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.epi_cfg = self.config.get("epidemiology", {})
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "epi"
        )

    def run(self, who_df: pd.DataFrame) -> dict:
        """Run full epidemiology analysis."""
        print("  [Epi] Analyzing TB trends...")
        trends = self.compute_trends(who_df)
        
        print("  [Epi] Detecting MDR patterns...")
        mdr = self.mdr_patterns(who_df)
        
        print("  [Epi] Aggregating district-level data...")
        district = self.district_aggregation(who_df)
        
        # Save results
        safe_save_csv(trends, self.output_dir / "tb_trends.csv")
        safe_save_csv(mdr, self.output_dir / "mdr_patterns.csv")
        safe_save_csv(district, self.output_dir / "district_aggregation.csv")
        
        summary = {
            "n_regions": who_df["Region"].nunique() if "Region" in who_df.columns else 0,
            "year_range": f"{who_df['Year'].min()}-{who_df['Year'].max()}" if "Year" in who_df.columns else "N/A",
            "overall_incidence_trend": self._trend_direction(trends, "TB_Incidence_per_100k"),
            "overall_mdr_trend": self._trend_direction(mdr, "MDR_TB_Percentage"),
            "highest_burden_region": district.iloc[0]["Region"] if len(district) > 0 else "N/A",
            "national_mean_mdr_pct": round(who_df["MDR_TB_Percentage"].mean(), 2) if "MDR_TB_Percentage" in who_df.columns else 0,
        }
        safe_save_json(summary, self.output_dir / "epi_summary.json")
        
        return {
            "trends": trends,
            "mdr_patterns": mdr,
            "district_aggregation": district,
            "summary": summary,
        }

    def compute_trends(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute temporal trends for each region using linear regression."""
        results = []
        
        for region in df["Region"].unique():
            region_data = df[df["Region"] == region].sort_values("Year")
            years = region_data["Year"].values
            
            for metric in ["TB_Incidence_per_100k", "TB_Mortality_per_100k", "MDR_TB_Percentage"]:
                if metric not in region_data.columns:
                    continue
                values = region_data[metric].values
                
                if len(years) >= 3:
                    slope, intercept, r_val, p_val, std_err = stats.linregress(years, values)
                    # Annual percent change
                    mean_val = np.mean(values)
                    apc = (slope / mean_val * 100) if mean_val != 0 else 0
                else:
                    slope = intercept = r_val = p_val = std_err = apc = 0
                
                results.append({
                    "Region": region,
                    "Metric": metric,
                    "Slope": round(slope, 4),
                    "Intercept": round(intercept, 2),
                    "R_squared": round(r_val ** 2, 4),
                    "P_value": round(p_val, 6),
                    "Annual_Pct_Change": round(apc, 2),
                    "Trend": "Increasing" if slope > 0 else "Decreasing",
                    "Start_Value": round(values[0], 2),
                    "End_Value": round(values[-1], 2),
                })
        
        result_df = pd.DataFrame(results)
        print(f"  [Epi] ✓ Computed trends for {df['Region'].nunique()} regions × {len(result_df) // max(df['Region'].nunique(), 1)} metrics")
        return result_df

    def mdr_patterns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Analyze MDR-TB prevalence patterns and risk classification."""
        results = []
        
        for region in df["Region"].unique():
            region_data = df[df["Region"] == region].sort_values("Year")
            
            if "MDR_TB_Percentage" not in region_data.columns:
                continue
            
            mdr_vals = region_data["MDR_TB_Percentage"].values
            years = region_data["Year"].values
            
            mean_mdr = np.mean(mdr_vals)
            latest_mdr = mdr_vals[-1] if len(mdr_vals) > 0 else 0
            max_mdr = np.max(mdr_vals)
            
            # Trend analysis
            if len(years) >= 3:
                slope, _, _, p_val, _ = stats.linregress(years, mdr_vals)
            else:
                slope, p_val = 0, 1.0
            
            # Risk classification
            if latest_mdr > 10:
                risk = "Critical"
            elif latest_mdr > 6:
                risk = "High"
            elif latest_mdr > 3:
                risk = "Moderate"
            else:
                risk = "Low"
            
            # Volatility (coefficient of variation)
            cv = (np.std(mdr_vals) / mean_mdr * 100) if mean_mdr > 0 else 0
            
            results.append({
                "Region": region,
                "Mean_MDR_Pct": round(mean_mdr, 2),
                "Latest_MDR_Pct": round(latest_mdr, 2),
                "Max_MDR_Pct": round(max_mdr, 2),
                "MDR_Trend_Slope": round(slope, 4),
                "Trend_P_value": round(p_val, 6),
                "Trend_Direction": "Rising" if slope > 0.1 else ("Falling" if slope < -0.1 else "Stable"),
                "Risk_Category": risk,
                "Volatility_CV": round(cv, 2),
            })
        
        result_df = pd.DataFrame(results)
        result_df = result_df.sort_values("Latest_MDR_Pct", ascending=False).reset_index(drop=True)
        
        n_critical = (result_df["Risk_Category"] == "Critical").sum()
        n_high = (result_df["Risk_Category"] == "High").sum()
        print(f"  [Epi] ✓ MDR patterns: {n_critical} critical, {n_high} high-risk regions")
        return result_df

    def district_aggregation(self, df: pd.DataFrame) -> pd.DataFrame:
        """Aggregate and rank regions by TB burden composite score."""
        agg = df.groupby("Region").agg({
            "TB_Incidence_per_100k": ["mean", "std", "max"],
            "TB_Mortality_per_100k": ["mean", "max"],
            "MDR_TB_Percentage": ["mean", "max"],
            "Cases_Notified": ["sum", "mean"] if "Cases_Notified" in df.columns else ["count"],
            "Treatment_Success_Rate": ["mean"] if "Treatment_Success_Rate" in df.columns else ["count"],
        })
        
        # Flatten column names
        agg.columns = ["_".join(col).strip() for col in agg.columns.values]
        agg = agg.reset_index()
        
        # Composite burden score (normalized weighted sum)
        for col in agg.columns[1:]:
            max_val = agg[col].max()
            if max_val > 0:
                agg[f"{col}_norm"] = agg[col] / max_val
        
        norm_cols = [c for c in agg.columns if c.endswith("_norm")]
        if norm_cols:
            agg["Burden_Score"] = agg[norm_cols].mean(axis=1)
        else:
            agg["Burden_Score"] = 0
            
        agg = agg.sort_values("Burden_Score", ascending=False).reset_index(drop=True)
        agg["Burden_Rank"] = range(1, len(agg) + 1)
        
        print(f"  [Epi] ✓ District aggregation: {len(agg)} regions ranked")
        return agg

    def _trend_direction(self, df: pd.DataFrame, metric: str) -> str:
        """Get overall trend direction for a metric."""
        if "Metric" in df.columns:
            subset = df[df["Metric"] == metric]
            if len(subset) > 0:
                mean_slope = subset["Slope"].mean()
                return "Increasing" if mean_slope > 0 else "Decreasing"
        elif "MDR_Trend_Slope" in df.columns:
            mean_slope = df["MDR_Trend_Slope"].mean()
            return "Rising" if mean_slope > 0 else "Falling"
        return "Stable"


if __name__ == "__main__":
    from utils.helpers import generate_mock_who_data
    engine = EpidemiologyEngine()
    who_df = generate_mock_who_data()
    results = engine.run(who_df)
    print(f"\nEpi Summary: {results['summary']}")
