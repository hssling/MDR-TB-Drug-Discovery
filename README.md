---
title: MDR-TB Drug Discovery Pipeline
emoji: 🧬
colorFrom: indigo
colorTo: red
sdk: gradio
app_file: challenges/huggingface/app.py
pinned: false
---

# MDR-TB AI Drug Discovery Pipeline

An end-to-end computational platform for multidrug-resistant tuberculosis drug discovery, resistance-aware target prioritization, generative chemistry, structural screening, and automated research-document generation.

This project combines data ingestion, target scoring, docking, post-docking triage, machine learning, manuscript generation, Hugging Face challenge deployment, and Kaggle-oriented dataset packaging in one reproducible workflow.

---

## Why This Project Matters

MDR-TB remains a major therapeutic challenge. This repository is designed to help computationally:

- prioritize resistance-relevant TB targets;
- generate and rank novel candidate compounds;
- triage leads using docking, ADMET, MD proxy, off-target, quantum, and retrosynthesis modules;
- generate manuscript and proposal outputs from saved artefacts; and
- share the workflow through Hugging Face and Kaggle.

---

## Core Capabilities

### Multi-Stage Discovery Pipeline

- **Resistance-aware scoring** in `models/resistance_module.py`
- **Target prioritization** in `engines/target_engine.py`
- **De novo compound generation** in `models/de_novo_generator.py`
- **Docking and structural screening** in `engines/docking_engine.py`
- **Post-docking triage** in `engines/md_engine.py`, `models/admet_predictor.py`, `engines/polypharmacology_engine.py`, `engines/qmmm_engine.py`, and `engines/retrosynthesis_engine.py`
- **Machine-learning support** in `models/ml_pipeline.py` and `models/ranking_engine.py`

### Automated Research Outputs

- **Publication-ready manuscript generation** in `generators/manuscript_generator.py`
- **ICMR proposal generation** in `generators/icmr_generator.py`
- **IEC protocol generation** in `generators/iec_generator.py`

### Community and Deployment Hooks

- **Hugging Face challenge prototype** in `challenges/huggingface/app.py`
- **Kaggle challenge and starter assets** in `challenges/kaggle/`
- **Dashboards** in `dashboards/`

---

## Repository Map

| Area | Purpose |
| :--- | :--- |
| `run_pipeline.py` | Main orchestration entry point |
| `data_connectors/` | Data loaders and external-source connectors |
| `engines/` | Core scientific analysis modules |
| `models/` | Predictive, ranking, and generation models |
| `generators/` | Manuscript, ICMR, and IEC document generation |
| `dashboards/` | Interactive dashboards |
| `challenges/huggingface/` | Hugging Face Space app |
| `challenges/kaggle/` | Kaggle challenge documentation and starter script |
| `outputs/` | Generated artefacts and reports |
| `scripts/prepare_release_bundle.py` | Clean export builder for public releases |

---

## Quick Start

### Local Run

```bash
pip install -r requirements.txt
python run_pipeline.py
```

Outputs are written under `outputs/`.

### Docker

```bash
docker-compose up --build
```

This launches the pipeline plus two Streamlit dashboards exposed on `8501` and `8502` via `docker-compose.yml`.

---

## Output Scope

The current workflow can generate:

- ranked targets and ranked compounds;
- docking reports and post-docking summaries;
- ADMET, MD proxy, off-target, QM, and retrosynthesis outputs;
- manuscript-ready markdown files;
- ICMR and IEC draft documents; and
- summary JSON files for downstream sharing.

Example generated artefacts include:

- `outputs/manuscript/manuscript_v7_final.md`
- `outputs/docking/docking_results_InhA.csv`
- `outputs/md_simulations/md_summary.csv`
- `outputs/admet/admet_toxicity_report.csv`
- `outputs/polypharmacology/human_offtarget_report.csv`
- `outputs/quantum_mechanics/electronic_orbitals.csv`
- `outputs/retrosynthesis/synthesis_routes.csv`

---

## Hugging Face Scope

This repository is configured to support a Hugging Face Space using the metadata block at the top of this `README.md`.

### Included Space Capability

- **Interactive SMILES challenge app** via `challenges/huggingface/app.py`
- **Gradio-based prototype evaluation flow** for user-submitted compounds
- **Public-facing demonstration layer** for sharing the project beyond the codebase

### Intended Hugging Face Use

- demo the project interactively;
- let users explore the compound challenge interface;
- showcase the broader computational drug-discovery narrative; and
- provide a lightweight public engagement front end.

### Current Space Limitations

- the Space app is a challenge/demo layer, not the full pipeline UI;
- it evaluates user SMILES with a lightweight ChemBERTa-based flow rather than rerunning the full docking pipeline; and
- it should be interpreted as an interactive outreach component, not a scientific validation interface.

---

## Kaggle Scope

This repository also supports Kaggle dataset publishing and challenge-oriented sharing.

### Included Kaggle Assets

- `dataset-metadata.json`
- `challenges/kaggle/KAGGLE_COMPETITION.md`
- `challenges/kaggle/starter_script.py`

### Intended Kaggle Use

- publish cleaned project snapshots as datasets;
- distribute generated outputs for reproducible analysis;
- host notebooks that benchmark or extend the ranking workflow; and
- support competition-style community experiments around candidate prioritization.

### Clean Kaggle Publishing Workflow

To avoid uploading `.git`, `__pycache__`, `.pyc`, and other local-only artefacts:

```bash
python scripts/prepare_release_bundle.py --output-dir release_bundle
kaggle datasets version -p release_bundle --dir-mode zip -m "Update dataset"
```

### Current CI Behavior

The repository CI workflow in `.github/workflows/ci.yml` currently deploys Kaggle from the repository root with `--dir-mode zip`. The clean-bundle script above is the recommended manual publishing path when you want a tidier public dataset.

---

## Challenges and Competitions

This repository includes public-facing challenge components intended to stimulate community exploration.

### Hugging Face Challenge

Located in `challenges/huggingface/app.py`, this interactive app lets users submit candidate SMILES strings and compare them conceptually against the internal pipeline lead through a lightweight challenge interface.

### Kaggle Challenge Assets

Located in `challenges/kaggle/`, these files provide:

- challenge framing;
- dataset-oriented onboarding;
- starter code for notebook-based exploration; and
- a path toward reproducible community benchmarking.

These challenge components are for computational exploration and engagement, not for making validated efficacy claims.

---

## Publication and Documentation

The repository can generate and maintain publication-oriented outputs directly from the computational workflow. This includes:

- full manuscript drafts;
- supplementary materials;
- proposal drafts; and
- ethics/exemption text.

This makes the repository useful not only as a scientific codebase, but also as a documentation and dissemination pipeline.

---

## Safety and Scope Boundaries

This project is for computational discovery support only.

It does **not** provide:

- wet-lab protocols;
- operational synthesis instructions;
- biological assay procedures;
- clinical guidance; or
- validated therapeutic claims.

All computational findings require appropriate expert review and experimental validation.

---

## Release Notes for Public Publishing

For public pushes:

- Hugging Face uses the metadata in this `README.md`
- Kaggle exports should be built from `scripts/prepare_release_bundle.py`
- generated outputs in `outputs/` can be shared selectively depending on the release target

---

---

## Drug Discovery Agent (AI-Powered, Disease-Agnostic)

A fully autonomous, Claude-orchestrated agent for end-to-end computational drug discovery — for **any disease**, not just TB. Give it a natural-language prompt; it searches real databases, computes ADMET properties, trains QSAR models, ranks compounds, and writes a full report.

### Why This Agent?

All results come exclusively from real, published data (ChEMBL, PubChem, PDB, EuropePMC/PubMed). Built-in guardrails reject fabricated data, missing citations, and low-confidence results at every step.

### Quick Start

```bash
# Install agent dependencies
pip install -r requirements_agent.txt

# Set your Anthropic API key
export ANTHROPIC_API_KEY=sk-ant-...   # or copy .env.example → .env

# Run a drug discovery task
python -m agent "Find novel DprE1 inhibitors for drug-resistant tuberculosis"
python -m agent "Identify EGFR kinase inhibitors for lung cancer"
python -m agent "Discover ACE2 inhibitors for COVID-19"
python -m agent "Find PfDHFR inhibitors for malaria"
```

Reports are saved to `outputs/agent_runs/<run_id>/` as `report.md`, `report.json`, and `compounds.csv`.

### Available Tools (12)

| Tool | Source | Description |
|------|--------|-------------|
| `search_chembl_target` | ChEMBL | Find protein targets by name/organism |
| `fetch_chembl_activities` | ChEMBL | Retrieve IC50/Ki/Kd bioactivity data |
| `get_chembl_compound` | ChEMBL | Full molecule record (SMILES, drug status) |
| `search_pubchem_compound` | PubChem | Search compounds by name or CID |
| `get_pubchem_bioassays` | PubChem | BioAssay activity summary for a CID |
| `search_pdb_structure` | RCSB PDB | Find protein structures by name |
| `download_pdb_structure` | RCSB PDB | Download PDB file + extract metadata |
| `search_literature` | EuropePMC | Search peer-reviewed publications |
| `fetch_paper_details` | EuropePMC | Full abstract + metadata by PMID/DOI |
| `compute_admet` | RDKit (local) | MW, cLogP, TPSA, QED, Lipinski, hERG |
| `train_qsar_model` | RDKit + scikit-learn | Binary classifier on ECFP4 fingerprints |
| `rank_compounds` | RDKit (local) | Composite score: pIC50 + LBVS + QED + Lipinski |

### MCP Server (Claude Desktop / Cursor Integration)

```bash
# stdio transport (Claude Desktop)
python -m agent.mcp

# HTTP transport (web MCP clients)
python -m agent.mcp --http --port 8765
```

**Claude Desktop config** (`claude_desktop_config.json`):
```json
{
  "mcpServers": {
    "drug-discovery": {
      "command": "python",
      "args": ["-m", "agent.mcp"],
      "cwd": "<path to this repo>"
    }
  }
}
```

### CLI Options

```
python -m agent "prompt" [--model MODEL] [--max-iterations N] [--output-dir DIR] [--quiet]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--model` | `claude-sonnet-4-6` | Claude model ID |
| `--max-iterations` | 25 | Agentic loop cap |
| `--output-dir` | `outputs/agent_runs` | Report output directory |
| `--quiet` | off | Suppress per-tool progress |

### Python API

```python
from agent import DrugDiscoveryAgent, AgentConfig

cfg = AgentConfig.from_env()
agent = DrugDiscoveryAgent(config=cfg)
report = agent.run("Find InhA inhibitors for MDR-TB")
report.print_summary()
# Access: report.top_compounds, report.report_path, report.metadata
```

### Guardrails

Every tool result is validated before being returned to Claude:

- **Hard failures**: missing citations on success, confidence outside [0, 1], fabrication patterns (random/mock data), success=False with no error message
- **Soft warnings**: confidence < 0.30, small training datasets (< 50 compounds)

---

## Contact

**Dr Siddalingaiah H S**
Professor, Community Medicine
Shridevi Institute of Medical Sciences and Research Hospital, Tumkur
Email: `hssling@yahoo.com`
ORCID: `0000-0002-4771-8285`
