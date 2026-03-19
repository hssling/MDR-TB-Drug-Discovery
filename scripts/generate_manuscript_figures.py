from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"
FIGURES = OUTPUTS / "figures"


def setup() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_workflow() -> None:
    steps = [
        "Data loading",
        "Omics + epidemiology",
        "Resistance mapping",
        "Target ranking",
        "De novo generation",
        "Docking + ranking",
        "ADMET / MD / QM / polypharm",
        "Retrosynthesis + manuscript",
    ]
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis("off")
    x_positions = [0.02, 0.14, 0.27, 0.39, 0.52, 0.64, 0.77, 0.90]
    for x, label in zip(x_positions, steps):
        box = FancyBboxPatch((x - 0.055, 0.35), 0.11, 0.28, boxstyle="round,pad=0.02", facecolor="#dbeafe", edgecolor="#1d4ed8", linewidth=1.5)
        ax.add_patch(box)
        ax.text(x, 0.49, label, ha="center", va="center", fontsize=10, wrap=True)
    for start, end in zip(x_positions[:-1], x_positions[1:]):
        ax.annotate("", xy=(end - 0.07, 0.49), xytext=(start + 0.07, 0.49), arrowprops=dict(arrowstyle="->", lw=1.8, color="#1f2937"))
    ax.set_title("Figure 1. Computational workflow for MDR-TB lead prioritization", fontsize=14, weight="bold")
    save(fig, "figure_1_workflow.png")


def figure_target_ranking() -> None:
    df = pd.read_csv(OUTPUTS / "targets" / "scored_targets.csv").sort_values("Final_Score", ascending=True).tail(10)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(df["Target"], df["Final_Score"], color="#2563eb")
    ax.set_xlabel("Final score")
    ax.set_title("Figure 2. Top-ranked TB drug targets")
    save(fig, "figure_2_target_ranking.png")


def figure_mdr_patterns() -> None:
    df = pd.read_csv(OUTPUTS / "epi" / "mdr_patterns.csv").sort_values("Latest_MDR_Pct", ascending=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(df["Region"], df["Latest_MDR_Pct"], color="#dc2626")
    ax.set_ylabel("Latest MDR-TB %")
    ax.set_title("Figure 3. Regional MDR-TB burden patterns")
    ax.tick_params(axis="x", rotation=45)
    save(fig, "figure_3_mdr_patterns.png")


def figure_pathways() -> None:
    df = pd.read_csv(OUTPUTS / "omics" / "pathway_scores.csv").head(8).sort_values("Enrichment_Score", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(df["Pathway"], df["Enrichment_Score"], color="#059669")
    ax.set_xlabel("Enrichment score")
    ax.set_title("Figure 4. Top enriched pathways in the omics module")
    save(fig, "figure_4_pathways.png")


def figure_docking() -> None:
    df = pd.read_csv(OUTPUTS / "docking" / "docking_results_InhA.csv").head(10).sort_values("Binding_Affinity_kcal_mol", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(df["Compound_ID"], df["Binding_Affinity_kcal_mol"], color="#7c3aed")
    ax.set_xlabel("Binding affinity (kcal/mol)")
    ax.set_title("Figure 5. Top docking hits against InhA")
    save(fig, "figure_5_docking_hits.png")


def figure_molecular_structures() -> None:
    ranked = pd.read_csv(OUTPUTS / "ranking" / "ranked_compounds.csv")
    top = ranked[ranked["Compound_ID"].notna()][["Compound_ID", "SMILES"]]
    mols = []
    legends = []
    for _, row in top.iterrows():
        mol = Chem.MolFromSmiles(row["SMILES"])
        if mol is not None:
            mols.append(mol)
            legends.append(row["Compound_ID"])
        if len(mols) == 6:
            break
    image = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300, 250), legends=legends, useSVG=False)
    image.save(FIGURES / "figure_6_top_structures.png")


def figure_lead_profile() -> None:
    admet = pd.read_csv(OUTPUTS / "admet" / "admet_toxicity_report.csv")
    md = pd.read_csv(OUTPUTS / "md_simulations" / "md_summary.csv")
    qm = pd.read_csv(OUTPUTS / "quantum_mechanics" / "electronic_orbitals.csv")
    poly = pd.read_csv(OUTPUTS / "polypharmacology" / "human_offtarget_report.csv")

    lead = admet[admet["Compound_ID"] == "MDR_AI_030"].iloc[0]
    md_row = md.iloc[0]
    qm_row = qm.iloc[0]
    poly_row = poly.iloc[0]

    labels = ["QED", "1-Tanimoto", "RMSD score", "H-bonds/5", "Band gap/5"]
    values = [
        float(lead["QED_Drug_Likeness"]),
        1 - float(poly_row["Max_Tanimoto_to_Human_Toxin"]),
        max(0, 1 - float(md_row["RMSD_nm"])),
        min(1, float(md_row["Persistent_H_Bonds"]) / 5.0),
        min(1, float(qm_row["Band_Gap_eV"]) / 5.0),
    ]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(labels, values, color=["#0ea5e9", "#10b981", "#f59e0b", "#8b5cf6", "#ef4444"])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Scaled lead-profile score")
    ax.set_title("Figure 7. Integrated computational profile of MDR_AI_030")
    ax.tick_params(axis="x", rotation=20)
    save(fig, "figure_7_lead_profile.png")


def main() -> None:
    setup()
    figure_workflow()
    figure_target_ranking()
    figure_mdr_patterns()
    figure_pathways()
    figure_docking()
    figure_molecular_structures()
    figure_lead_profile()
    print(f"Generated figures in: {FIGURES}")


if __name__ == "__main__":
    main()
