from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import Draw


ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "outputs"
FIGURES = OUTPUTS / "figures"


def setup() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    RDLogger.DisableLog("rdApp.*")


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def annotate_horizontal_bars(ax: plt.Axes, values: list[float], fmt: str = "{:.2f}", pad: float = 0.02) -> None:
    xmin, xmax = ax.get_xlim()
    span = xmax - xmin
    for patch, value in zip(ax.patches, values):
        y = patch.get_y() + patch.get_height() / 2
        x = patch.get_width()
        if value >= 0:
            text_x = x + span * pad
            ha = "left"
        else:
            text_x = x - span * pad
            ha = "right"
        ax.text(text_x, y, fmt.format(value), va="center", ha=ha, fontsize=9, color="#111827")


def annotate_vertical_bars(ax: plt.Axes, values: list[float], fmt: str = "{:.2f}", pad: float = 0.02) -> None:
    ymin, ymax = ax.get_ylim()
    span = ymax - ymin
    for patch, value in zip(ax.patches, values):
        x = patch.get_x() + patch.get_width() / 2
        y = patch.get_height()
        ax.text(x, y + span * pad, fmt.format(value), va="bottom", ha="center", fontsize=9, color="#111827")


def figure_workflow() -> None:
    steps = [
        "Data loading",
        "Omics + epidemiology",
        "Resistance mapping",
        "Target ranking",
        "De novo chemistry",
        "Docking + ranking",
        "ADMET / MD / QM /\noff-target screening",
        "Retrosynthesis +\nmanuscript outputs",
    ]
    fig, ax = plt.subplots(figsize=(8.5, 12))
    ax.axis("off")
    positions = [(0.5, 0.9 - index * 0.105) for index in range(len(steps))]
    for index, ((x, y), label) in enumerate(zip(positions, steps), start=1):
        box = FancyBboxPatch(
            (x - 0.28, y - 0.038),
            0.56,
            0.076,
            boxstyle="round,pad=0.02",
            facecolor="#dbeafe",
            edgecolor="#1d4ed8",
            linewidth=1.5,
        )
        ax.add_patch(box)
        ax.text(x - 0.24, y, f"{index}", ha="center", va="center", fontsize=10, weight="bold", color="#1e3a8a")
        ax.text(x - 0.2, y, label, ha="left", va="center", fontsize=11, weight="bold")

    for start, end in zip(positions[:-1], positions[1:]):
        ax.annotate(
            "",
            xy=(end[0], end[1] + 0.05),
            xytext=(start[0], start[1] - 0.05),
            arrowprops=dict(arrowstyle="->", lw=1.8, color="#1f2937"),
        )

    ax.set_title("Figure 1. Computational workflow for MDR-TB lead prioritization", fontsize=15, weight="bold", pad=14)
    save(fig, "figure_1_workflow.png")


def figure_target_ranking() -> None:
    df = pd.read_csv(OUTPUTS / "targets" / "scored_targets.csv").sort_values("Final_Score", ascending=False).head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    values = df["Final_Score"].tolist()
    ax.barh(df["Target"], values, color="#2563eb")
    ax.set_xlabel("Final score")
    ax.set_title("Figure 2. Top-ranked TB drug targets")
    ax.set_xlim(0, max(values) + 0.06)
    annotate_horizontal_bars(ax, values)
    save(fig, "figure_2_target_ranking.png")


def figure_mdr_patterns() -> None:
    df = pd.read_csv(OUTPUTS / "epi" / "mdr_patterns.csv").sort_values("Latest_MDR_Pct", ascending=False).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    values = df["Latest_MDR_Pct"].tolist()
    ax.barh(df["Region"], values, color="#dc2626")
    ax.set_xlabel("Latest MDR-TB %")
    ax.set_title("Figure 3. Regional MDR-TB burden patterns")
    ax.set_xlim(0, max(values) + 1.0)
    annotate_horizontal_bars(ax, values, fmt="{:.1f}")
    save(fig, "figure_3_mdr_patterns.png")


def figure_pathways() -> None:
    df = pd.read_csv(OUTPUTS / "omics" / "pathway_scores.csv").head(8).sort_values("Enrichment_Score", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 5.6))
    values = df["Enrichment_Score"].tolist()
    ax.barh(df["Pathway"], values, color="#059669")
    ax.set_xlabel("Enrichment score")
    ax.set_title("Figure 4. Top enriched pathways in the omics module")
    ax.set_xlim(0, max(values) + 0.6)
    annotate_horizontal_bars(ax, values)
    save(fig, "figure_4_pathways.png")


def figure_docking() -> None:
    df = pd.read_csv(OUTPUTS / "docking" / "docking_results_InhA.csv").sort_values("Binding_Affinity_kcal_mol", ascending=True).head(10)
    fig, ax = plt.subplots(figsize=(9, 5.6))
    values = df["Binding_Affinity_kcal_mol"].tolist()
    ax.barh(df["Compound_ID"], values, color="#7c3aed")
    ax.set_xlabel("Binding affinity (kcal/mol)")
    ax.set_title("Figure 5. Top docking hits against InhA")
    ax.set_xlim(min(values) - 0.5, 0)
    annotate_horizontal_bars(ax, values)
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

    labels = ["QED", "1-Tanimoto", "RMSD\nscore", "H-bonds/5", "Band gap/5"]
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
    annotate_vertical_bars(ax, values)
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
