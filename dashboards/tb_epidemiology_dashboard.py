"""
TB Epidemiology Dashboard — Phase 10
======================================
Interactive Streamlit dashboard for visualizing
TB incidence, mortality, and MDR patterns across Indian states.

SAFETY: Computational only — visualization of public health data.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys
import json

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

st.set_page_config(
    page_title="TB Epidemiology Dashboard — MDR-TB AI Pipeline v3",
    page_icon="🦠",
    layout="wide",
    initial_sidebar_state="expanded"
)


def load_data():
    """Load epidemiology data from pipeline outputs."""
    base = Path(__file__).parent.parent
    
    # Load WHO data
    who_path = base / "data" / "who" / "who_tb_data.csv"
    if who_path.exists():
        who_df = pd.read_csv(who_path)
    else:
        # Generate mock data
        from utils.helpers import generate_mock_who_data
        who_df = generate_mock_who_data()
    
    # Load analysis results
    trends_path = base / "outputs" / "epi" / "tb_trends.csv"
    trends_df = pd.read_csv(trends_path) if trends_path.exists() else None
    
    mdr_path = base / "outputs" / "epi" / "mdr_patterns.csv"
    mdr_df = pd.read_csv(mdr_path) if mdr_path.exists() else None
    
    district_path = base / "outputs" / "epi" / "district_aggregation.csv"
    district_df = pd.read_csv(district_path) if district_path.exists() else None
    
    summary_path = base / "outputs" / "epi" / "epi_summary.json"
    summary = {}
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)
    
    return who_df, trends_df, mdr_df, district_df, summary


def main():
    # Custom CSS
    st.markdown("""
    <style>
    .main-header { font-size: 2.5rem; color: #1e3a5f; font-weight: 700; margin-bottom: 0; }
    .sub-header { font-size: 1.1rem; color: #666; margin-top: 0; }
    .metric-card { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                   color: white; padding: 20px; border-radius: 12px; text-align: center; }
    .warning-card { background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
                    color: white; padding: 20px; border-radius: 12px; }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<p class="main-header">🦠 TB Epidemiology Dashboard</p>', unsafe_allow_html=True)
    st.markdown('<p class="sub-header">MDR-TB AI Pipeline v3 — Public Health Surveillance Analytics</p>', unsafe_allow_html=True)
    st.markdown("---")
    
    who_df, trends_df, mdr_df, district_df, summary = load_data()
    
    # Sidebar filters
    st.sidebar.title("🔍 Filters")
    regions = st.sidebar.multiselect(
        "Select Regions",
        options=sorted(who_df["Region"].unique()),
        default=sorted(who_df["Region"].unique())[:5]
    )
    
    year_range = st.sidebar.slider(
        "Year Range",
        min_value=int(who_df["Year"].min()),
        max_value=int(who_df["Year"].max()),
        value=(int(who_df["Year"].min()), int(who_df["Year"].max()))
    )
    
    # Filter data
    filtered = who_df[
        (who_df["Region"].isin(regions)) &
        (who_df["Year"] >= year_range[0]) &
        (who_df["Year"] <= year_range[1])
    ]
    
    # Key metrics
    st.subheader("📊 Key Indicators")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Mean Incidence", f"{filtered['TB_Incidence_per_100k'].mean():.1f}/100k",
                   delta=f"{summary.get('overall_incidence_trend', 'N/A')}")
    with col2:
        st.metric("Mean Mortality", f"{filtered['TB_Mortality_per_100k'].mean():.1f}/100k")
    with col3:
        st.metric("Mean MDR %", f"{filtered['MDR_TB_Percentage'].mean():.1f}%")
    with col4:
        total_cases = filtered["Cases_Notified"].sum() if "Cases_Notified" in filtered.columns else 0
        st.metric("Total Cases Notified", f"{total_cases:,.0f}")
    
    st.markdown("---")
    
    # Charts
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Incidence Trends", "🔴 MDR Patterns", "🗺️ Regional Burden", "📋 Data Table"])
    
    with tab1:
        st.subheader("TB Incidence Trends Over Time")
        fig1 = px.line(
            filtered, x="Year", y="TB_Incidence_per_100k",
            color="Region", markers=True,
            title="TB Incidence per 100,000 Population",
            template="plotly_white"
        )
        fig1.update_layout(height=500, xaxis_title="Year", yaxis_title="Incidence per 100k")
        st.plotly_chart(fig1, use_container_width=True)
        
        # Mortality trends
        fig2 = px.line(
            filtered, x="Year", y="TB_Mortality_per_100k",
            color="Region", markers=True,
            title="TB Mortality per 100,000 Population",
            template="plotly_white"
        )
        fig2.update_layout(height=400)
        st.plotly_chart(fig2, use_container_width=True)
    
    with tab2:
        st.subheader("MDR-TB Prevalence Patterns")
        fig3 = px.line(
            filtered, x="Year", y="MDR_TB_Percentage",
            color="Region", markers=True,
            title="MDR-TB Percentage by Region",
            template="plotly_white"
        )
        fig3.update_layout(height=500)
        st.plotly_chart(fig3, use_container_width=True)
        
        # MDR Risk Heatmap
        if mdr_df is not None:
            st.subheader("MDR Risk Classification")
            risk_colors = {"Critical": "🔴", "High": "🟠", "Moderate": "🟡", "Low": "🟢"}
            mdr_display = mdr_df.copy()
            if "Risk_Category" in mdr_display.columns:
                mdr_display["Risk_Icon"] = mdr_display["Risk_Category"].map(risk_colors)
                st.dataframe(mdr_display[["Region", "Risk_Icon", "Risk_Category", "Latest_MDR_Pct", "Trend_Direction"]],
                            use_container_width=True)
    
    with tab3:
        st.subheader("Regional Burden Analysis")
        
        # Bar chart: latest year data
        latest_year = filtered["Year"].max()
        latest_data = filtered[filtered["Year"] == latest_year]
        
        fig4 = px.bar(
            latest_data.sort_values("TB_Incidence_per_100k", ascending=True),
            x="TB_Incidence_per_100k", y="Region",
            orientation="h",
            color="MDR_TB_Percentage",
            color_continuous_scale="Reds",
            title=f"TB Incidence & MDR Burden by Region ({latest_year})",
            template="plotly_white"
        )
        fig4.update_layout(height=500)
        st.plotly_chart(fig4, use_container_width=True)
        
        # Scatter: Incidence vs MDR
        fig5 = px.scatter(
            latest_data, x="TB_Incidence_per_100k", y="MDR_TB_Percentage",
            size="Cases_Notified" if "Cases_Notified" in latest_data.columns else None,
            color="Region",
            title="Incidence vs MDR Percentage",
            template="plotly_white"
        )
        fig5.update_layout(height=400)
        st.plotly_chart(fig5, use_container_width=True)
    
    with tab4:
        st.subheader("Raw Data")
        st.dataframe(filtered, use_container_width=True)
        st.download_button(
            "📥 Download Filtered Data",
            filtered.to_csv(index=False),
            "tb_epi_filtered.csv",
            "text/csv"
        )
    
    # Footer
    st.markdown("---")
    st.caption("MDR-TB AI Pipeline v3 — Computational analysis only. Data sources: WHO Global TB Programme (mock/synthetic for demonstration).")


if __name__ == "__main__":
    main()
