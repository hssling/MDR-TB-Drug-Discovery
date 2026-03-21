---
title: MDR-TB Drug Discovery Pipeline
emoji: 🧬
colorFrom: indigo
colorTo: red
sdk: gradio
app_file: app.py
pinned: false
---

# MDR-TB Drug Discovery — AI Pipeline + Drug Discovery Agent

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![ChEMBL](https://img.shields.io/badge/data-ChEMBL-green.svg)](https://www.ebi.ac.uk/chembl/)
[![Kaggle Dataset](https://img.shields.io/badge/dataset-Kaggle-20BEFF.svg)](https://www.kaggle.com/datasets/jkhospital/mdrtb-drug-discovery-ai-pipeline)
[![HuggingFace Space](https://img.shields.io/badge/demo-HuggingFace-yellow.svg)](https://huggingface.co/spaces/hssling/MDR-TB-Drug-Discovery)

> **Two tools in one repo:**
> 1. **MDR-TB Computational Pipeline** — end-to-end resistance-aware drug discovery for tuberculosis
> 2. **Disease-Agnostic AI Drug Discovery Agent** — Claude-orchestrated, real-data-only agent for *any* disease

---

## Key Results (MDR-TB Pipeline)

| Metric | Value |
|--------|-------|
| ChEMBL target | CHEMBL1849 (InhA, *M. tuberculosis*) |
| Compounds screened | 421 with measured IC50 |
| Top-ranked compound | CHEMBL3125270 — IC50 4 nM, Composite Score 0.825 |
| QSAR AUC (5-fold CV) | 0.979 (Random Forest) |
| Lipinski-compliant leads | 100% of top 10 |
| katG S315T resistance rate | ~67% (CRyPTIC n=12,289) |
| WHO burden (2023) | 7.8 M cases, 1.13 M deaths |

All compound scores use **real, measured ChEMBL IC50 values** — no computational docking predictions are used as ground truth.

---

## What's in This Repository

```
.
├── agent/                    ← AI Drug Discovery Agent (any disease)
├── challenges/
│   ├── huggingface/app.py    ← Gradio demo (InhA challenge)
│   └── kaggle/               ← Competition docs + starter code
├── engines/                  ← Core scientific modules (docking, MD, ADMET…)
├── generators/               ← Manuscript, proposal, ethics doc generation
├── models/                   ← QSAR, ranking, de novo generation
├── outputs/
│   ├── ranking/ranked_compounds.csv   ← Primary dataset (421 compounds)
│   ├── manuscript/manuscript_v8_genuine.md
│   └── figures/              ← Publication figures
├── scripts/                  ← Export, submission, release bundle helpers
├── requirements.txt          ← Pipeline dependencies
├── requirements_agent.txt    ← Agent dependencies
└── .env.example              ← API key template
```

---

## Part 1 — Drug Discovery Agent (Any Disease)

An autonomous Claude-orchestrated agent that conducts full drug discovery workflows from a single natural-language prompt. Designed for researchers studying any disease target.

### Install

```bash
pip install -r requirements_agent.txt
cp .env.example .env        # add your ANTHROPIC_API_KEY
```

Get an API key at [console.anthropic.com](https://console.anthropic.com).

### Run

```bash
# Tuberculosis
python -m agent "Find novel DprE1 inhibitors for drug-resistant tuberculosis"

# Other diseases
python -m agent "Identify EGFR kinase inhibitors for lung cancer"
python -m agent "Discover ACE2 inhibitors for COVID-19"
python -m agent "Find PfDHFR inhibitors for malaria"
python -m agent "Identify PCSK9 inhibitors for cardiovascular disease"
```

Reports are saved to `outputs/agent_runs/<run_id>/` as `report.md`, `report.json`, and `compounds.csv`.

### CLI Options

```
python -m agent "prompt" [options]

Options:
  --model MODEL           Claude model ID (default: claude-sonnet-4-6)
  --max-iterations N      Agentic loop cap (default: 25)
  --output-dir DIR        Report output directory (default: outputs/agent_runs)
  --quiet, -q             Suppress per-tool progress output
```

### Python API

```python
from agent import DrugDiscoveryAgent, AgentConfig

cfg = AgentConfig.from_env()
agent = DrugDiscoveryAgent(config=cfg)
report = agent.run("Find InhA inhibitors for MDR-TB", verbose=True)
report.print_summary()

# Access structured results
print(report.top_compounds)     # list of ranked compound dicts
print(report.report_path)       # path to report.md
print(report.metadata)          # run metadata (model, iterations, timestamp)
```

### MCP Server (Claude Desktop / Cursor)

Exposes all 12 tools via the [Model Context Protocol](https://modelcontextprotocol.io) so any MCP-compatible AI client can use them.

```bash
# stdio transport — for Claude Desktop
python -m agent.mcp

# HTTP transport — for web MCP clients
python -m agent.mcp --http --port 8765
```

**Claude Desktop config** (`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS, `%APPDATA%\Claude\claude_desktop_config.json` on Windows):

```json
{
  "mcpServers": {
    "drug-discovery": {
      "command": "python",
      "args": ["-m", "agent.mcp"],
      "cwd": "/path/to/this/repo",
      "env": { "ANTHROPIC_API_KEY": "sk-ant-..." }
    }
  }
}
```

### Available Tools

| # | Tool | Data Source | Description |
|---|------|------------|-------------|
| 1 | `search_chembl_target` | ChEMBL | Find protein targets by name and organism |
| 2 | `fetch_chembl_activities` | ChEMBL | Retrieve IC50/Ki/Kd bioactivity data |
| 3 | `get_chembl_compound` | ChEMBL | Full molecule record: SMILES, drug status, synonyms |
| 4 | `search_pubchem_compound` | PubChem | Search compounds by name or CID |
| 5 | `get_pubchem_bioassays` | PubChem | BioAssay activity summary for a compound |
| 6 | `search_pdb_structure` | RCSB PDB | Find protein crystal structures |
| 7 | `download_pdb_structure` | RCSB PDB | Download PDB file + extract ligands and resolution |
| 8 | `search_literature` | EuropePMC | Search peer-reviewed publications |
| 9 | `fetch_paper_details` | EuropePMC | Full abstract + metadata by PMID or DOI |
| 10 | `compute_admet` | RDKit (local) | MW, cLogP, TPSA, QED, Lipinski, hERG flag |
| 11 | `train_qsar_model` | RDKit + scikit-learn | Binary classifier (ECFP4 fingerprints, 5-fold CV) |
| 12 | `rank_compounds` | RDKit (local) | Composite score: pIC50 (40%) + LBVS (30%) + QED (20%) + Lipinski (10%) |

### Guardrails

Every tool result passes through `GuardrailEnforcer` before being returned to Claude:

| Check | Type | Action on failure |
|-------|------|-------------------|
| Missing citations on success | Hard | Raise `GuardrailViolation` |
| Confidence outside [0, 1] | Hard | Raise `GuardrailViolation` |
| Fabrication patterns (random/mock data) | Hard | Raise `GuardrailViolation` |
| success=False with no error message | Hard | Raise `GuardrailViolation` |
| Confidence < 0.30 | Soft | Add warning, continue |
| Dataset < 50 compounds for QSAR | Soft | Add warning, continue |

---

## Part 2 — MDR-TB Computational Pipeline

A full multi-stage pipeline for tuberculosis drug discovery specifically:

### Pipeline Stages

1. **Data ingestion** — WHO, CRyPTIC, ChEMBL, GEO omics data (`data_connectors/`)
2. **Resistance scoring** — katG S315T, rpoB frequency-weighted target scoring (`models/resistance_module.py`)
3. **Target prioritization** — multi-criteria scoring of druggable TB targets (`engines/target_engine.py`)
4. **Ligand-based virtual screening** — Tanimoto + ECFP4 vs known InhA inhibitors (`engines/docking_engine.py`)
5. **ADMET filtering** — RDKit Lipinski, Veber, QED, hERG (`models/admet_predictor.py`)
6. **QSAR modeling** — Random Forest / Gradient Boosting on ECFP4 fingerprints (`models/ml_pipeline.py`)
7. **Composite ranking** — 0.40×pIC50 + 0.30×LBVS + 0.20×QED + 0.10×Lipinski (`models/ranking_engine.py`)
8. **Document generation** — manuscript, ICMR proposal, IEC protocol (`generators/`)

### Quick Start

```bash
pip install -r requirements.txt
python run_pipeline.py
```

Outputs are written to `outputs/`.

### Docker

```bash
docker-compose up --build
```

Launches the pipeline plus two Streamlit dashboards on ports 8501 and 8502.

---

## Primary Dataset

**`outputs/ranking/ranked_compounds.csv`** — 421 InhA inhibitors from ChEMBL1849

Key columns:

| Column | Description |
|--------|-------------|
| `molecule_chembl_id` | ChEMBL compound identifier |
| `canonical_smiles` | RDKit-canonical SMILES string |
| `IC50_nM` | Measured IC50 in nM (from ChEMBL — not predicted) |
| `pIC50` | −log10(IC50 in M) |
| `MW`, `LogP`, `TPSA`, `HBD`, `HBA` | RDKit physicochemical properties |
| `QED` | Drug-likeness score [0–1] (Bickerton 2012) |
| `Lipinski_Pass` | Lipinski Rule-of-Five compliance |
| `hERG_Risk` | hERG cardiotoxicity flag (Waring 2010) |
| `LBVS_Composite_Score` | Tanimoto + pharmacophore similarity to known inhibitors |
| `QSAR_Active_Prob` | Random Forest predicted activity probability |
| `Composite_Score` | Final ranking score [0–1] |
| `Rank` | Overall rank (1 = best) |

---

## Community Resources

| Resource | Link |
|----------|------|
| GitHub | [hssling/MDR-TB-Drug-Discovery](https://github.com/hssling/MDR-TB-Drug-Discovery) |
| HuggingFace Space | [hssling/MDR-TB-Drug-Discovery](https://huggingface.co/spaces/hssling/MDR-TB-Drug-Discovery) |
| Kaggle Dataset | [jkhospital/mdrtb-drug-discovery-ai-pipeline](https://www.kaggle.com/datasets/jkhospital/mdrtb-drug-discovery-ai-pipeline) |
| Manuscript | Submitted to *Tuberculosis* (Elsevier, ISSN 1472-9792) |

---

## How to Contribute

1. Fork this repository
2. Run the pipeline or agent on a new target
3. Open an issue or PR with your findings
4. Share your Kaggle notebook — see [challenges/kaggle/](challenges/kaggle/)

---

## Citation

If you use this pipeline or dataset in your research, please cite:

```bibtex
@article{siddalingaiah2025mdrtb,
  title   = {Ligand-Based Virtual Screening of InhA Inhibitors for
             Multidrug-Resistant Tuberculosis Using ChEMBL Bioactivity Data},
  author  = {Siddalingaiah, H S},
  journal = {Tuberculosis},
  year    = {2025},
  note    = {Submitted},
  url     = {https://github.com/hssling/MDR-TB-Drug-Discovery}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

All ChEMBL data is used under the [ChEMBL Terms of Use](https://www.ebi.ac.uk/chembl/about/).
PDB structures are used under the [RCSB PDB Terms](https://www.rcsb.org/pages/usage-policy).

---

## Contact

**Dr Siddalingaiah H S**
Professor, Community Medicine
Shridevi Institute of Medical Sciences and Research Hospital, Tumkur, Karnataka, India
Email: `hssling@yahoo.com`
ORCID: [0000-0002-4771-8285](https://orcid.org/0000-0002-4771-8285)

---

## Safety Notice

This project is for **computational research support only**. It does not provide:

- wet-lab protocols or synthesis instructions
- clinical guidance or therapeutic recommendations
- validated efficacy or safety claims

All computational findings require appropriate expert review and experimental validation before any use in drug development.
