"""
Omics Engine — Phase 3
=======================
Differential expression analysis, pathway scoring, and gene ranking
for MDR-TB transcriptomic data.

SAFETY: Computational only — statistical analysis of expression data.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.helpers import load_config, ensure_dir, safe_save_csv, safe_save_json


class OmicsEngine:
    """Perform differential expression, pathway scoring, and gene ranking."""

    # TB-relevant KEGG pathways with associated gene indices (mock mapping)
    TB_PATHWAYS = {
        "mtu00010_Glycolysis": list(range(0, 30)),
        "mtu00020_TCA_cycle": list(range(30, 55)),
        "mtu00190_Oxidative_phosphorylation": list(range(55, 90)),
        "mtu00230_Purine_metabolism": list(range(90, 130)),
        "mtu00240_Pyrimidine_metabolism": list(range(130, 160)),
        "mtu00550_Peptidoglycan_biosynthesis": list(range(160, 185)),
        "mtu00620_Pyruvate_metabolism": list(range(185, 210)),
        "mtu00640_Propanoate_metabolism": list(range(210, 230)),
        "mtu00071_Fatty_acid_degradation": list(range(230, 260)),
        "mtu00361_Mycolic_acid_biosynthesis": list(range(260, 300)),
        "mtu02020_Two_component_system": list(range(300, 340)),
        "mtu03010_Ribosome": list(range(340, 380)),
        "mtu03030_DNA_replication": list(range(380, 400)),
    }

    def __init__(self, config: dict = None):
        self.config = config or load_config()
        self.omics_cfg = self.config.get("omics", {})
        self.log2fc_thresh = self.omics_cfg.get("log2fc_threshold", 1.0)
        self.pval_thresh = self.omics_cfg.get("pvalue_threshold", 0.05)
        self.top_n = self.omics_cfg.get("top_n_genes", 50)
        self.output_dir = ensure_dir(
            Path(self.config["paths"]["output_dir"]) / "omics"
        )

    def run(self, expression_df: pd.DataFrame) -> dict:
        """Run full omics analysis pipeline."""
        print("  [Omics] Running differential expression analysis...")
        de_results = self.differential_expression(expression_df)
        
        print("  [Omics] Scoring pathways...")
        pathway_scores = self.pathway_scoring(de_results, expression_df)
        
        print("  [Omics] Ranking genes...")
        gene_ranks = self.gene_ranking(de_results)
        
        # Save results
        safe_save_csv(de_results, self.output_dir / "differential_expression.csv")
        safe_save_csv(pathway_scores, self.output_dir / "pathway_scores.csv")
        safe_save_csv(gene_ranks, self.output_dir / "gene_rankings.csv")
        
        summary = {
            "total_genes_tested": len(de_results),
            "significant_up": int(((de_results["log2FC"] > self.log2fc_thresh) & 
                                  (de_results["pvalue"] < self.pval_thresh)).sum()),
            "significant_down": int(((de_results["log2FC"] < -self.log2fc_thresh) & 
                                    (de_results["pvalue"] < self.pval_thresh)).sum()),
            "top_pathway": pathway_scores.iloc[0]["Pathway"] if len(pathway_scores) > 0 else "N/A",
            "top_gene": gene_ranks.iloc[0]["Gene"] if len(gene_ranks) > 0 else "N/A",
        }
        safe_save_json(summary, self.output_dir / "omics_summary.json")
        
        return {
            "de_results": de_results,
            "pathway_scores": pathway_scores,
            "gene_rankings": gene_ranks,
            "summary": summary,
        }

    def differential_expression(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """
        Perform differential expression analysis (MDR vs Control).
        Uses Welch's t-test with Benjamini-Hochberg correction.
        """
        ctrl_cols = [c for c in expr_df.columns if c.startswith("CTRL")]
        mdr_cols = [c for c in expr_df.columns if c.startswith("MDR")]
        
        if not ctrl_cols or not mdr_cols:
            # Fallback: split columns in half
            mid = len(expr_df.columns) // 2
            ctrl_cols = expr_df.columns[:mid].tolist()
            mdr_cols = expr_df.columns[mid:].tolist()
        
        results = []
        for gene in expr_df.index:
            ctrl_vals = expr_df.loc[gene, ctrl_cols].values.astype(float)
            mdr_vals = expr_df.loc[gene, mdr_cols].values.astype(float)
            
            mean_ctrl = np.mean(ctrl_vals)
            mean_mdr = np.mean(mdr_vals)
            log2fc = mean_mdr - mean_ctrl  # Already in log2 space
            
            # Welch's t-test
            t_stat, pval = stats.ttest_ind(mdr_vals, ctrl_vals, equal_var=False)
            if np.isnan(pval):
                pval = 1.0
            
            results.append({
                "Gene": gene,
                "mean_CTRL": round(mean_ctrl, 4),
                "mean_MDR": round(mean_mdr, 4),
                "log2FC": round(log2fc, 4),
                "t_statistic": round(t_stat, 4),
                "pvalue": pval,
            })
        
        df = pd.DataFrame(results)
        
        # Benjamini-Hochberg FDR correction
        df = df.sort_values("pvalue")
        n = len(df)
        df["rank"] = range(1, n + 1)
        df["padj"] = df["pvalue"] * n / df["rank"]
        df["padj"] = df["padj"].clip(upper=1.0)
        # Ensure monotonicity
        df["padj"] = df["padj"][::-1].cummin()[::-1]
        
        df["significant"] = (df["padj"] < self.pval_thresh) & (abs(df["log2FC"]) > self.log2fc_thresh)
        df = df.sort_values("padj")
        
        n_sig = df["significant"].sum()
        print(f"  [Omics] ✓ DE analysis: {n_sig} significant genes (padj < {self.pval_thresh}, |log2FC| > {self.log2fc_thresh})")
        
        return df

    def pathway_scoring(self, de_results: pd.DataFrame, expr_df: pd.DataFrame) -> pd.DataFrame:
        """
        Score TB-relevant pathways using gene set enrichment approach.
        Uses mean absolute log2FC of pathway genes as activity score.
        """
        pathway_data = []
        gene_list = de_results["Gene"].tolist()
        
        for pathway_name, gene_indices in self.TB_PATHWAYS.items():
            # Map indices to actual genes
            pathway_genes = [gene_list[i] for i in gene_indices if i < len(gene_list)]
            pathway_de = de_results[de_results["Gene"].isin(pathway_genes)]
            
            if len(pathway_de) == 0:
                continue
            
            mean_fc = pathway_de["log2FC"].abs().mean()
            n_sig = pathway_de["significant"].sum()
            mean_pval = pathway_de["padj"].mean()
            
            # Combined pathway score
            enrichment_score = mean_fc * (1 + n_sig / max(len(pathway_de), 1))
            
            # Permutation-based p-value estimate
            n_perm = 1000
            np.random.seed(hash(pathway_name) % 2**31)
            perm_scores = []
            for _ in range(n_perm):
                perm_genes = de_results.sample(n=len(pathway_de))
                perm_score = perm_genes["log2FC"].abs().mean() * (1 + perm_genes["significant"].sum() / max(len(perm_genes), 1))
                perm_scores.append(perm_score)
            
            pathway_pval = (np.sum(np.array(perm_scores) >= enrichment_score) + 1) / (n_perm + 1)
            
            pathway_data.append({
                "Pathway": pathway_name,
                "N_Genes": len(pathway_de),
                "N_Significant": int(n_sig),
                "Mean_abs_log2FC": round(mean_fc, 4),
                "Enrichment_Score": round(enrichment_score, 4),
                "Pathway_Pvalue": round(pathway_pval, 4),
            })
        
        df = pd.DataFrame(pathway_data)
        df = df.sort_values("Enrichment_Score", ascending=False).reset_index(drop=True)
        print(f"  [Omics] ✓ Pathway scoring: {len(df)} pathways scored")
        return df

    def gene_ranking(self, de_results: pd.DataFrame) -> pd.DataFrame:
        """
        Rank genes by a composite score combining fold change,
        statistical significance, and effect size.
        """
        df = de_results.copy()
        
        # Composite ranking score
        df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
        df["abs_log2FC"] = df["log2FC"].abs()
        
        # Normalize both to [0, 1]
        max_nlp = df["neg_log10_padj"].max()
        max_fc = df["abs_log2FC"].max()
        
        if max_nlp > 0:
            df["norm_significance"] = df["neg_log10_padj"] / max_nlp
        else:
            df["norm_significance"] = 0
            
        if max_fc > 0:
            df["norm_effect"] = df["abs_log2FC"] / max_fc
        else:
            df["norm_effect"] = 0
        
        # Composite score (weighted harmonic mean)
        w_sig, w_eff = 0.5, 0.5
        df["composite_score"] = (
            (w_sig + w_eff) /
            (w_sig / (df["norm_significance"] + 1e-10) + w_eff / (df["norm_effect"] + 1e-10))
        )
        
        df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
        df["rank"] = range(1, len(df) + 1)
        
        top_genes = df.head(self.top_n)[["rank", "Gene", "log2FC", "padj", "composite_score"]]
        print(f"  [Omics] ✓ Gene ranking: Top gene = {df.iloc[0]['Gene']} (score={df.iloc[0]['composite_score']:.4f})")
        
        return df[["rank", "Gene", "log2FC", "padj", "norm_significance", "norm_effect", "composite_score"]]


if __name__ == "__main__":
    from utils.helpers import generate_mock_expression_matrix
    engine = OmicsEngine()
    expr = generate_mock_expression_matrix()
    results = engine.run(expr)
    print(f"\nOmics Summary: {results['summary']}")
