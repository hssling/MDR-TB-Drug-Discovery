"""
Drug Discovery Dashboard — Phase 10
=====================================
Interactive Streamlit dashboard for visualizing compound analysis,
ML model results, target scores, and ranking results.

SAFETY: Computational only — visualization of in-silico results.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys
import json

sys.path.insert(0, str(Path(__file__).parent.parent))

st.set_page_config(
    page_title="Drug Discovery Dashboard — MDR-TB AI Pipeline v3",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_results():
    """Load pipeline results."""
    base = Path(__file__).parent.parent
    data = {}
    
    file_map = {
        "compounds": "outputs/compounds/descriptors.csv",
        "clustered": "outputs/compounds/clustered_compounds.csv",
        "lipinski": "outputs/compounds/lipinski_compliance.csv",
        "targets": "outputs/targets/scored_targets.csv",
        "rankings": "outputs/ranking/ranked_compounds.csv",
        "top10": "outputs/ranking/top_10_compounds.csv",
        "models": "outputs/models/model_comparison.csv",
        "features": "outputs/models/feature_importance.csv",
        "resistance": "outputs/resistance/resistance_scores.csv",
        "mutations": "outputs/resistance/mutation_catalog.csv",
        "de_results": "outputs/omics/differential_expression.csv",
        "gene_ranks": "outputs/omics/gene_rankings.csv",
        "pathways": "outputs/omics/pathway_scores.csv",
    }
    
    for key, path in file_map.items():
        full_path = base / path
        if full_path.exists():
            data[key] = pd.read_csv(full_path)
        else:
            data[key] = None
    
    # Load summaries
    for summary_name in ["compound_summary", "target_summary", "ranking_summary", "ml_summary"]:
        summary_path = base / "outputs" / summary_name.replace("_summary", "") / f"{summary_name}.json"
        # Try alternative paths
        for alt_dir in ["compounds", "targets", "ranking", "models"]:
            alt_path = base / "outputs" / alt_dir / f"{summary_name}.json"
            if alt_path.exists():
                with open(alt_path) as f:
                    data[summary_name] = json.load(f)
                break
        else:
            data[summary_name] = {}
    
    return data


def main():
    st.markdown("""
    <style>
    .main-header { font-size: 2.5rem; color: #1e3a5f; font-weight: 700; }
    .sub-header { font-size: 1.1rem; color: #666; }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<p class="main-header">💊 Drug Discovery Dashboard</p>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">MDR-TB AI Pipeline v3 — Computational Drug Discovery Analytics</p>', unsafe_allow_html=True)
    st.markdown("---")
    
    data = load_results()
    
    # Navigation tabs
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "🧬 Omics", "🎯 Targets", "💊 Compounds", "🤖 ML Models", "🛡️ Resistance", "🏆 Final Ranking"
    ])
    
    # ---- TAB 1: OMICS ----
    with tab1:
        st.subheader("🧬 Differential Expression & Pathway Analysis")
        
        if data["de_results"] is not None:
            de = data["de_results"]
            col1, col2, col3 = st.columns(3)
            with col1:
                n_sig = de["significant"].sum() if "significant" in de.columns else 0
                st.metric("Significant Genes", n_sig)
            with col2:
                n_up = ((de["log2FC"] > 1) & (de.get("significant", False))).sum() if "log2FC" in de.columns else 0
                st.metric("Upregulated", n_up)
            with col3:
                n_down = ((de["log2FC"] < -1) & (de.get("significant", False))).sum() if "log2FC" in de.columns else 0
                st.metric("Downregulated", n_down)
            
            # Volcano plot
            if "log2FC" in de.columns and "padj" in de.columns:
                de_plot = de.copy()
                de_plot["neg_log10_padj"] = -np.log10(de_plot["padj"].clip(lower=1e-300))
                de_plot["Status"] = "Not Significant"
                de_plot.loc[(de_plot["log2FC"] > 1) & (de_plot["padj"] < 0.05), "Status"] = "Upregulated"
                de_plot.loc[(de_plot["log2FC"] < -1) & (de_plot["padj"] < 0.05), "Status"] = "Downregulated"
                
                fig_volcano = px.scatter(
                    de_plot.head(1000), x="log2FC", y="neg_log10_padj",
                    color="Status",
                    color_discrete_map={"Upregulated": "#e74c3c", "Downregulated": "#3498db", "Not Significant": "#95a5a6"},
                    title="Volcano Plot — Differential Expression",
                    template="plotly_white", opacity=0.6
                )
                fig_volcano.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
                fig_volcano.add_vline(x=1, line_dash="dash", line_color="gray")
                fig_volcano.add_vline(x=-1, line_dash="dash", line_color="gray")
                fig_volcano.update_layout(height=500)
                st.plotly_chart(fig_volcano, use_container_width=True)
        
        if data["pathways"] is not None:
            st.subheader("Pathway Enrichment Scores")
            fig_path = px.bar(
                data["pathways"].head(10), x="Enrichment_Score", y="Pathway",
                orientation="h", color="Enrichment_Score",
                color_continuous_scale="Viridis",
                title="Top 10 Enriched Pathways",
                template="plotly_white"
            )
            fig_path.update_layout(height=400)
            st.plotly_chart(fig_path, use_container_width=True)
    
    # ---- TAB 2: TARGETS ----
    with tab2:
        st.subheader("🎯 Drug Target Scoring")
        if data["targets"] is not None:
            targets = data["targets"]
            
            # Radar chart for top targets
            top_targets = targets.head(5)
            fig_radar = go.Figure()
            categories = ["Druggability", "Essentiality", "Conservation", "Resistance_Association", "Expression_FC_Score"]
            available_cats = [c for c in categories if c in top_targets.columns]
            
            for _, row in top_targets.iterrows():
                values = [row[c] for c in available_cats]
                values.append(values[0])  # Close the polygon
                fig_radar.add_trace(go.Scatterpolar(
                    r=values,
                    theta=available_cats + [available_cats[0]],
                    fill="toself",
                    name=row["Target"],
                    opacity=0.6
                ))
            
            fig_radar.update_layout(
                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                title="Top 5 Drug Targets — Multi-Factor Profile",
                height=500
            )
            st.plotly_chart(fig_radar, use_container_width=True)
            
            # Target bar chart
            fig_targets = px.bar(
                targets, x="Final_Score", y="Target",
                color="Priority",
                color_discrete_map={"High": "#e74c3c", "Medium": "#f39c12", "Low": "#95a5a6"},
                orientation="h",
                title="All Targets Ranked by Score",
                template="plotly_white"
            )
            fig_targets.update_layout(height=500)
            st.plotly_chart(fig_targets, use_container_width=True)
    
    # ---- TAB 3: COMPOUNDS ----
    with tab3:
        st.subheader("💊 Compound Analysis")
        if data["clustered"] is not None:
            compounds = data["clustered"]
            
            col1, col2 = st.columns(2)
            with col1:
                if "LogP" in compounds.columns and "MolWt" in compounds.columns:
                    fig_scatter = px.scatter(
                        compounds, x="MolWt", y="LogP",
                        color="Cluster" if "Cluster" in compounds.columns else None,
                        title="Chemical Space — MolWt vs LogP",
                        template="plotly_white"
                    )
                    fig_scatter.update_layout(height=400)
                    st.plotly_chart(fig_scatter, use_container_width=True)
            
            with col2:
                if "TPSA" in compounds.columns:
                    fig_hist = px.histogram(
                        compounds, x="TPSA", color="Cluster" if "Cluster" in compounds.columns else None,
                        title="TPSA Distribution by Cluster",
                        template="plotly_white"
                    )
                    fig_hist.update_layout(height=400)
                    st.plotly_chart(fig_hist, use_container_width=True)
            
            # Lipinski compliance
            if data["lipinski"] is not None:
                lip = data["lipinski"]
                if "Lipinski_Pass" in lip.columns:
                    pass_count = lip["Lipinski_Pass"].sum()
                    fail_count = len(lip) - pass_count
                    fig_lip = px.pie(
                        values=[pass_count, fail_count],
                        names=["Pass", "Fail"],
                        title="Lipinski Rule of Five Compliance",
                        color_discrete_sequence=["#2ecc71", "#e74c3c"]
                    )
                    st.plotly_chart(fig_lip, use_container_width=True)
    
    # ---- TAB 4: ML MODELS ----
    with tab4:
        st.subheader("🤖 Machine Learning Model Results")
        if data["models"] is not None:
            models = data["models"]
            
            # Model comparison bar chart
            metrics_to_plot = ["accuracy", "precision", "recall", "f1_score", "roc_auc"]
            available_metrics = [m for m in metrics_to_plot if m in models.columns]
            
            if available_metrics and "algorithm" in models.columns:
                fig_models = go.Figure()
                for metric in available_metrics:
                    fig_models.add_trace(go.Bar(
                        name=metric.replace("_", " ").title(),
                        x=models["algorithm"],
                        y=models[metric]
                    ))
                fig_models.update_layout(
                    barmode="group",
                    title="Model Performance Comparison",
                    template="plotly_white",
                    height=400
                )
                st.plotly_chart(fig_models, use_container_width=True)
            
            st.dataframe(models, use_container_width=True)
        
        # Feature importance
        if data["features"] is not None:
            fig_fi = px.bar(
                data["features"], x="Importance", y="Feature",
                orientation="h", title="Feature Importance",
                template="plotly_white", color="Importance",
                color_continuous_scale="Blues"
            )
            st.plotly_chart(fig_fi, use_container_width=True)
    
    # ---- TAB 5: RESISTANCE ----
    with tab5:
        st.subheader("🛡️ MDR Resistance Analysis")
        if data["resistance"] is not None:
            resistance = data["resistance"]
            
            fig_res = px.bar(
                resistance, x="Resistance_Score", y="Gene",
                color="Priority",
                color_discrete_map={"Critical": "#e74c3c", "High": "#f39c12", "Moderate": "#3498db"},
                orientation="h",
                title="Resistance Gene Scores",
                template="plotly_white"
            )
            fig_res.update_layout(height=400)
            st.plotly_chart(fig_res, use_container_width=True)
        
        if data["mutations"] is not None:
            st.subheader("Mutation Catalog")
            st.dataframe(data["mutations"], use_container_width=True)
    
    # ---- TAB 6: FINAL RANKING ----
    with tab6:
        st.subheader("🏆 Final Compound Ranking")
        if data["rankings"] is not None:
            rankings = data["rankings"]
            
            # Top 10 showcase
            top = rankings.head(10)
            score_cols = [c for c in ["score_activity", "score_lipinski", "score_target", "score_novelty", "score_resistance"] if c in top.columns]
            
            if score_cols and "Compound_ID" in top.columns:
                fig_rank = go.Figure()
                for col in score_cols:
                    fig_rank.add_trace(go.Bar(
                        name=col.replace("score_", "").title(),
                        x=top["Compound_ID"],
                        y=top[col]
                    ))
                fig_rank.update_layout(
                    barmode="stack",
                    title="Top 10 Compounds — Score Breakdown",
                    template="plotly_white",
                    height=500
                )
                st.plotly_chart(fig_rank, use_container_width=True)
            
            # Full ranking table
            display_cols = ["Rank", "Compound_ID", "Final_Rank_Score", "Priority_Class"]
            available_display = [c for c in display_cols if c in rankings.columns]
            if available_display:
                st.dataframe(rankings[available_display].head(20), use_container_width=True)
            
            st.download_button(
                "📥 Download Full Rankings",
                rankings.to_csv(index=False),
                "compound_rankings.csv",
                "text/csv"
            )
    
    # Footer
    st.markdown("---")
    st.caption("MDR-TB AI Pipeline v3 — Computational drug discovery. No wet-lab protocols. All results are in-silico predictions.")


if __name__ == "__main__":
    main()
