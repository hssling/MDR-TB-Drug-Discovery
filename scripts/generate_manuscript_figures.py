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


# ---------------------------------------------------------------------------
# Figure 1 — Innovative 5-phase horizontal pipeline (complete redesign)
# ---------------------------------------------------------------------------

def figure_workflow() -> None:
    """
    Innovative 5-phase horizontal pipeline figure for MDR-TB lead discovery.

    Layout
    ------
    • Five color-coded phase columns (left → right)
    • Phase header at top of each column
    • Four module boxes per column with name + metric subtitle
    • Horizontal arrows between columns at mid-height
    • Bottom summary strip showing the compound-narrowing funnel
    • Dashed separator between module area and summary strip
    """

    FIG_W, FIG_H = 18.0, 10.5
    MARGIN_L = 0.5
    PHASE_W = 2.8
    PHASE_GAP = 0.7
    COL_STEP = PHASE_W + PHASE_GAP   # 3.5

    # Vertical geometry
    TITLE_Y = FIG_H - 0.20
    HDR_TOP = FIG_H - 0.50
    HDR_H = 1.30
    HDR_BOT = HDR_TOP - HDR_H        # 8.70

    MOD_FIRST_TOP = HDR_BOT - 0.18   # 8.52
    MOD_H = 1.46
    MOD_GAP = 0.11
    N_MODULES = 4

    # bottom of last module
    MOD_LAST_BOT = MOD_FIRST_TOP - N_MODULES * MOD_H - (N_MODULES - 1) * MOD_GAP
    # = 8.52 - 5.84 - 0.33 = 2.35

    STRIP_H = 0.82
    STRIP_TOP = MOD_LAST_BOT - 0.28  # 2.07
    STRIP_BOT = STRIP_TOP - STRIP_H  # 1.25

    ARROW_Y = (MOD_FIRST_TOP + MOD_LAST_BOT) / 2  # mid of module area ≈ 5.44

    # Phase definitions
    phases = [
        dict(
            label="PHASE I",
            title="Disease\nContext",
            hdr_bg="#1e3a8a",
            mod_bg="#eff6ff",
            edge="#2563eb",
            txt="#1e3a8a",
            modules=[
                ("Data Loading", "15 targets · 55 compounds"),
                ("Omics Analysis", "5,000 genes tested"),
                ("Pathway Scoring", "Glycolysis top-ranked"),
                ("Epidemiology", "10 Indian regions · 2015–2023"),
            ],
        ),
        dict(
            label="PHASE II",
            title="Resistance &\nTarget Intelligence",
            hdr_bg="#7f1d1d",
            mod_bg="#fff1f2",
            edge="#dc2626",
            txt="#7f1d1d",
            modules=[
                ("Resistance Mapping", "10 genes · 40 mutations"),
                ("CRyPTIC Integration", "rpoB, katG top-scored"),
                ("AlphaFold Structures", "3D protein models built"),
                ("Target Scoring", "InhA  #1  (score 0.8695)"),
            ],
        ),
        dict(
            label="PHASE III",
            title="Chemical\nDiscovery",
            hdr_bg="#14532d",
            mod_bg="#f0fdf4",
            edge="#16a34a",
            txt="#14532d",
            modules=[
                ("De Novo Generation", "55 novel structures"),
                ("Lipinski Screening", "Drug-likeness filtered"),
                ("Docking vs InhA", "MDR_AI_030: -9.77 kcal/mol"),
                ("ML Ranking", "Random Forest  AUC 0.727"),
            ],
        ),
        dict(
            label="PHASE IV",
            title="Lead\nTriage",
            hdr_bg="#581c87",
            mod_bg="#faf5ff",
            edge="#7c3aed",
            txt="#581c87",
            modules=[
                ("MD Proxy (10 ns)", "RMSD 0.12 nm · 4 H-bonds"),
                ("ADMET Profiling", "GI: Yes  BBB: No  hERG: Low"),
                ("Polypharmacology", "Tanimoto 0.146 (selective)"),
                ("Quantum Chemistry", "Band gap 2.9 eV"),
            ],
        ),
        dict(
            label="PHASE V",
            title="Translation",
            hdr_bg="#78350f",
            mod_bg="#fffbeb",
            edge="#d97706",
            txt="#78350f",
            modules=[
                ("Retrosynthesis", "Sonogashira · 1 step"),
                ("Synthetic Feasibility", "Fragments commercially available"),
                ("Lead Nomination", "[#1]  MDR_AI_030"),
                ("Manuscript Output", "Audited · Reproducible"),
            ],
        ),
    ]

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.axis("off")
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, FIG_H)

    # ── Subtle page background ──
    bg = FancyBboxPatch(
        (0.08, 0.08), FIG_W - 0.16, FIG_H - 0.16,
        boxstyle="round,pad=0.05",
        facecolor="#f8fafc",
        edgecolor="#e2e8f0",
        linewidth=1.0,
        zorder=0,
    )
    ax.add_patch(bg)

    # ── Overall title ──
    TOTAL_W = len(phases) * PHASE_W + (len(phases) - 1) * PHASE_GAP
    ax.text(
        MARGIN_L + TOTAL_W / 2, TITLE_Y,
        "Figure 1.  Resistance-Informed Computational Pipeline for MDR-TB Lead Discovery",
        ha="center", va="top", fontsize=13, weight="bold", color="#0f172a",
    )

    # ── Phase columns ──
    for pi, phase in enumerate(phases):
        x_left = MARGIN_L + pi * COL_STEP
        x_right = x_left + PHASE_W
        x_mid = x_left + PHASE_W / 2

        # Phase header box
        hdr = FancyBboxPatch(
            (x_left, HDR_BOT), PHASE_W, HDR_H,
            boxstyle="round,pad=0.06",
            facecolor=phase["hdr_bg"],
            edgecolor="white",
            linewidth=1.8,
            zorder=3,
        )
        ax.add_patch(hdr)

        # Phase label (small, near top of header)
        ax.text(
            x_mid, HDR_TOP - 0.22,
            phase["label"],
            ha="center", va="center",
            fontsize=8, weight="bold",
            color="white", alpha=0.80,
            zorder=4,
        )
        # Phase title (larger, lower in header)
        ax.text(
            x_mid, HDR_BOT + 0.50,
            phase["title"],
            ha="center", va="center",
            fontsize=10.5, weight="bold",
            color="white", linespacing=1.3,
            zorder=4,
        )

        # Module boxes
        for mi, (mod_name, mod_desc) in enumerate(phase["modules"]):
            m_top = MOD_FIRST_TOP - mi * (MOD_H + MOD_GAP)
            m_bot = m_top - MOD_H
            pad_x = 0.08

            mod_box = FancyBboxPatch(
                (x_left + pad_x, m_bot),
                PHASE_W - 2 * pad_x,
                MOD_H,
                boxstyle="round,pad=0.05",
                facecolor=phase["mod_bg"],
                edgecolor=phase["edge"],
                linewidth=1.2,
                zorder=3,
            )
            ax.add_patch(mod_box)

            # Module name (bold, near top of box)
            ax.text(
                x_mid, m_top - 0.36,
                mod_name,
                ha="center", va="center",
                fontsize=9, weight="bold",
                color=phase["txt"],
                zorder=4,
            )
            # Module metric / description (italic, near bottom)
            ax.text(
                x_mid, m_bot + 0.38,
                mod_desc,
                ha="center", va="center",
                fontsize=7.5, style="italic",
                color="#374151",
                linespacing=1.25,
                zorder=4,
            )

        # Arrow to next phase
        if pi < len(phases) - 1:
            ax.annotate(
                "",
                xy=(x_right + PHASE_GAP - 0.10, ARROW_Y),
                xytext=(x_right + 0.10, ARROW_Y),
                arrowprops=dict(
                    arrowstyle="->,head_width=0.28,head_length=0.18",
                    lw=2.4,
                    color="#94a3b8",
                    connectionstyle="arc3,rad=0",
                ),
                zorder=5,
            )

    # ── Dashed separator ──
    sep_y = (MOD_LAST_BOT + STRIP_TOP) / 2
    ax.plot(
        [MARGIN_L, MARGIN_L + TOTAL_W],
        [sep_y, sep_y],
        color="#cbd5e1", lw=1.0, ls="--", zorder=2,
    )

    # ── Bottom summary strip (funnel) ──
    strip_items = [
        ("55 Compounds\nGenerated", "#2563eb"),
        ("15 Targets\nScored & Ranked", "#dc2626"),
        ("Top 10\nDocking Hits", "#16a34a"),
        ("Multi-module\nLead Triage", "#7c3aed"),
        ("[#1] MDR_AI_030\nNominated Lead", "#d97706"),
    ]

    # "Pipeline Summary" label on the left
    ax.text(
        MARGIN_L - 0.08,
        STRIP_BOT + STRIP_H / 2,
        "Pipeline\nSummary",
        ha="right", va="center",
        fontsize=7.5, weight="bold",
        color="#64748b", linespacing=1.3,
    )

    for si, (label, color) in enumerate(strip_items):
        sx = MARGIN_L + si * COL_STEP
        is_last = si == len(strip_items) - 1

        strip_box = FancyBboxPatch(
            (sx, STRIP_BOT), PHASE_W, STRIP_H,
            boxstyle="round,pad=0.05",
            facecolor=color,
            edgecolor="white",
            linewidth=1.5,
            alpha=0.92 if is_last else 0.55,
            zorder=3,
        )
        ax.add_patch(strip_box)

        ax.text(
            sx + PHASE_W / 2,
            STRIP_BOT + STRIP_H / 2,
            label,
            ha="center", va="center",
            fontsize=8.5,
            weight="bold" if is_last else "normal",
            color="white",
            linespacing=1.25,
            zorder=4,
        )

        # Arrow between strip items
        if si < len(strip_items) - 1:
            ax.annotate(
                "",
                xy=(sx + PHASE_W + PHASE_GAP - 0.10, STRIP_BOT + STRIP_H / 2),
                xytext=(sx + PHASE_W + 0.10, STRIP_BOT + STRIP_H / 2),
                arrowprops=dict(
                    arrowstyle="->,head_width=0.14,head_length=0.10",
                    lw=1.8,
                    color="#94a3b8",
                ),
                zorder=5,
            )

    fig.savefig(FIGURES / "figure_1_workflow.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 2 — Top enriched pathways (first plot in §3.1)
# ---------------------------------------------------------------------------

def figure_pathways() -> None:
    df = (
        pd.read_csv(OUTPUTS / "omics" / "pathway_scores.csv")
        .head(8)
        .sort_values("Enrichment_Score", ascending=True)
    )
    fig, ax = plt.subplots(figsize=(9, 5.6))
    values = df["Enrichment_Score"].tolist()
    ax.barh(df["Pathway"], values, color="#059669")
    ax.set_xlabel("Enrichment score")
    ax.set_title("Figure 2. Top enriched pathways in the omics module")
    ax.set_xlim(0, max(values) + 0.6)
    annotate_horizontal_bars(ax, values)
    save(fig, "figure_2_pathways.png")


# ---------------------------------------------------------------------------
# Figure 3 — Regional MDR-TB burden (second plot in §3.1)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Figure 4 — Target ranking (§3.2)
# ---------------------------------------------------------------------------

def figure_target_ranking() -> None:
    df = (
        pd.read_csv(OUTPUTS / "targets" / "scored_targets.csv")
        .sort_values("Final_Score", ascending=False)
        .head(10)
        .iloc[::-1]
    )
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    values = df["Final_Score"].tolist()
    ax.barh(df["Target"], values, color="#2563eb")
    ax.set_xlabel("Final score")
    ax.set_title("Figure 4. Top-ranked TB drug targets")
    ax.set_xlim(0, max(values) + 0.06)
    annotate_horizontal_bars(ax, values)
    save(fig, "figure_4_target_ranking.png")


# ---------------------------------------------------------------------------
# Figure 5 — Docking hits against InhA
# ---------------------------------------------------------------------------

def figure_docking() -> None:
    df = (
        pd.read_csv(OUTPUTS / "docking" / "docking_results_InhA.csv")
        .sort_values("Binding_Affinity_kcal_mol", ascending=True)
        .head(10)
    )
    fig, ax = plt.subplots(figsize=(9, 5.6))
    values = df["Binding_Affinity_kcal_mol"].tolist()
    ax.barh(df["Compound_ID"], values, color="#7c3aed")
    ax.set_xlabel("Binding affinity (kcal/mol)")
    ax.set_title("Figure 5. Top docking hits against InhA")
    ax.set_xlim(min(values) - 0.5, 0)
    annotate_horizontal_bars(ax, values)
    save(fig, "figure_5_docking_hits.png")


# ---------------------------------------------------------------------------
# Figure 6 — Top valid de novo molecular structures
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Figure 7 — Integrated lead profile for MDR_AI_030
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    setup()

    # Remove legacy filenames to avoid stale files alongside renamed ones
    for stale in ("figure_2_target_ranking.png", "figure_4_pathways.png"):
        p = FIGURES / stale
        if p.exists():
            p.unlink()

    figure_workflow()           # Figure 1 — innovative pipeline
    figure_pathways()           # Figure 2 — pathways (first in §3.1)
    figure_mdr_patterns()       # Figure 3 — MDR patterns (second in §3.1)
    figure_target_ranking()     # Figure 4 — target ranking (§3.2)
    figure_docking()            # Figure 5
    figure_molecular_structures()  # Figure 6
    figure_lead_profile()       # Figure 7

    print(f"Generated figures in: {FIGURES}")


if __name__ == "__main__":
    main()
