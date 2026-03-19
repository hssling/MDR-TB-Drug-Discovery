"""
MDR-TB AI Pipeline v3 — Utility Helpers
========================================
Common utility functions used across all pipeline modules.
SAFETY: Computational only — no wet-lab or synthesis steps.
"""

import os
import json
import yaml
import datetime
import pandas as pd
import numpy as np
from pathlib import Path


def load_config(config_path: str = None) -> dict:
    """Load YAML configuration file."""
    if config_path is None:
        config_path = Path(__file__).parent.parent / "config.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def ensure_dir(dir_path: str) -> Path:
    """Create directory if it doesn't exist and return Path."""
    p = Path(dir_path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def get_timestamp() -> str:
    """Return ISO-format timestamp string."""
    return datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def safe_save_csv(df: pd.DataFrame, path: str, **kwargs):
    """Save DataFrame to CSV, creating parent dirs if needed."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(p, index=False, **kwargs)
    print(f"  ✓ Saved: {p}")


def safe_save_json(data: dict, path: str):
    """Save dict to JSON, creating parent dirs if needed."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  ✓ Saved: {p}")


def log_phase(phase_name: str):
    """Print a formatted phase header."""
    bar = "=" * 60
    print(f"\n{bar}")
    print(f"  PHASE: {phase_name}")
    print(f"{bar}\n")


def generate_mock_expression_matrix(n_genes: int = 5000, n_samples: int = 20,
                                     seed: int = 42) -> pd.DataFrame:
    """Generate a synthetic gene expression matrix for testing."""
    np.random.seed(seed)
    genes = [f"MTB_gene_{i:04d}" for i in range(n_genes)]
    
    # Half control, half MDR-TB
    n_ctrl = n_samples // 2
    n_mdr = n_samples - n_ctrl
    samples_ctrl = [f"CTRL_{i+1}" for i in range(n_ctrl)]
    samples_mdr = [f"MDR_{i+1}" for i in range(n_mdr)]
    samples = samples_ctrl + samples_mdr
    
    # Base expression values (log2 scale)
    base = np.random.normal(loc=8, scale=2, size=(n_genes, n_samples))
    
    # Add differential expression for some genes (first 200 = DE genes)
    n_de = min(200, n_genes)
    fold_changes = np.random.choice([-1, 1], size=n_de) * np.random.uniform(1.5, 4.0, size=n_de)
    for i in range(n_de):
        base[i, n_ctrl:] += fold_changes[i]
    
    df = pd.DataFrame(base, index=genes, columns=samples)
    df.index.name = "Gene"
    return df


def generate_mock_who_data(regions: list = None, years: list = None,
                            seed: int = 42) -> pd.DataFrame:
    """Generate synthetic WHO TB incidence/mortality data."""
    np.random.seed(seed)
    if regions is None:
        regions = ["Maharashtra", "Gujarat", "Rajasthan", "UP", "TN",
                    "Karnataka", "Delhi", "WB", "MP", "Bihar"]
    if years is None:
        years = list(range(2015, 2024))
    
    rows = []
    for region in regions:
        base_incidence = np.random.randint(80, 300)
        base_mortality = np.random.randint(10, 60)
        base_mdr_pct = np.random.uniform(2.0, 12.0)
        for yr in years:
            trend = (yr - 2015) * np.random.uniform(-3, 1)
            rows.append({
                "Region": region,
                "Year": yr,
                "TB_Incidence_per_100k": max(10, base_incidence + trend + np.random.normal(0, 8)),
                "TB_Mortality_per_100k": max(1, base_mortality + trend * 0.3 + np.random.normal(0, 3)),
                "MDR_TB_Percentage": max(0.5, base_mdr_pct + (yr - 2015) * 0.3 + np.random.normal(0, 0.8)),
                "Treatment_Success_Rate": min(100, max(50, 75 + np.random.normal(0, 5))),
                "Cases_Notified": int(max(100, base_incidence * 50 + np.random.normal(0, 500))),
            })
    return pd.DataFrame(rows)


def generate_mock_compounds(n: int = 50, seed: int = 42) -> pd.DataFrame:
    """Generate synthetic compound data with SMILES-like identifiers."""
    np.random.seed(seed)
    # Real TB-relevant SMILES fragments + random extensions
    base_smiles = [
        "CC1=CC(=CC=C1)NC(=O)C2=CC=CC=C2",   # Generic amide
        "C1=CC=C(C=C1)C(=O)NN",                # Hydrazide
        "CC(=O)NNC(=S)NC1=CC=CC=C1",            # Thiosemicarbazone
        "C1=CC(=CC=C1)NC2=NC=NC3=CC=CC=C23",    # Quinazoline
        "CC1=CC2=C(C=C1)N=C(N2)SCC(=O)O",       # Benzimidazole
    ]
    
    rows = []
    for i in range(n):
        smiles = base_smiles[i % len(base_smiles)]
        rows.append({
            "Compound_ID": f"MDRTB_CPD_{i+1:04d}",
            "SMILES": smiles,
            "MolWt": np.random.uniform(150, 600),
            "LogP": np.random.uniform(-1, 5),
            "TPSA": np.random.uniform(20, 140),
            "NumHAcceptors": np.random.randint(1, 10),
            "NumHDonors": np.random.randint(0, 5),
            "NumRotatableBonds": np.random.randint(0, 10),
            "Activity_Label": np.random.choice([0, 1], p=[0.6, 0.4]),
            "Source": np.random.choice(["PubChem", "DrugBank", "Literature"]),
        })
    return pd.DataFrame(rows)


def generate_mock_drugbank_entries(n: int = 30, seed: int = 42) -> list:
    """Generate mock DrugBank-style entries for TB drugs."""
    np.random.seed(seed)
    tb_drugs = [
        ("Isoniazid", "DB00951", "C4H7N3O", "InhA"),
        ("Rifampicin", "DB01045", "C43H58N4O12", "RpoB"),
        ("Pyrazinamide", "DB00339", "C5H5N3O", "PncA"),
        ("Ethambutol", "DB00330", "C10H24N2O2", "EmbB"),
        ("Streptomycin", "DB01082", "C21H39N7O12", "rrs"),
        ("Levofloxacin", "DB01137", "C18H20FN3O4", "GyrA"),
        ("Moxifloxacin", "DB00218", "C21H24FN3O4", "GyrA"),
        ("Bedaquiline", "DB12207", "C32H31BrN2O2", "AtpE"),
        ("Delamanid", "DB12230", "C25H25F3N4O6", "F420"),
        ("Linezolid", "DB00601", "C16H20FN3O4", "rrl"),
    ]
    entries = []
    for i in range(min(n, len(tb_drugs) * 3)):
        base = tb_drugs[i % len(tb_drugs)]
        entries.append({
            "drugbank_id": base[1] if i < len(tb_drugs) else f"DB{90000+i:05d}",
            "name": base[0] if i < len(tb_drugs) else f"Compound_{i}",
            "formula": base[2],
            "target_gene": base[3],
            "mechanism": f"Inhibits {base[3]} in M. tuberculosis",
            "status": np.random.choice(["approved", "investigational", "experimental"]),
            "half_life": f"{np.random.uniform(1, 24):.1f} hours",
            "toxicity_class": np.random.choice(["low", "moderate", "high"]),
        })
    return entries
