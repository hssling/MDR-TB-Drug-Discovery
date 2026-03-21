# Drug Discovery Agent — Full Documentation

The Drug Discovery Agent is a Claude-orchestrated AI system for end-to-end computational drug discovery. It works for **any disease target** from a single natural-language prompt.

---

## Table of Contents

1. [Architecture](#architecture)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Tools Reference](#tools-reference)
5. [Guardrails](#guardrails)
6. [MCP Server](#mcp-server)
7. [Report Format](#report-format)
8. [Configuration](#configuration)
9. [Extending the Agent](#extending-the-agent)
10. [Troubleshooting](#troubleshooting)

---

## Architecture

```
User prompt
    │
    ▼
DrugDiscoveryOrchestrator
    │   Uses Anthropic API (claude-sonnet-4-6)
    │   Agentic tool_use loop (up to 25 iterations)
    │
    ├── ChEMBLTools     → ChEMBL REST API  (target search, bioactivity)
    ├── PubChemTools    → PubChem PUG REST (compound search, bioassays)
    ├── PDBTools        → RCSB PDB API     (structure search + download)
    ├── PubMedTools     → EuropePMC API    (literature search)
    └── RDKitTools      → Local computation (ADMET, QSAR, ranking)
            │
            ▼
    GuardrailEnforcer  ← validates every tool result before returning to Claude
            │
            ▼
    ReportGenerator    → report.md + report.json + compounds.csv
```

### Design Principles

- **Real data only** — all tools query live external databases or compute locally. No mock data, no hallucinated IC50 values.
- **Citations mandatory** — every successful tool result must include at least one `Citation` (source URL, database, year). The guardrail rejects results without citations.
- **Confidence-scored** — every result carries a float confidence in [0, 1] with a rationale string explaining the score.
- **Stateless tools** — each tool call is independent and idempotent. PDB files are cached locally.

---

## Installation

### Requirements

- Python 3.10+
- Anthropic API key ([console.anthropic.com](https://console.anthropic.com))
- Internet access (for ChEMBL, PubChem, PDB, EuropePMC)

### Steps

```bash
# Clone the repository
git clone https://github.com/hssling/MDR-TB-Drug-Discovery.git
cd MDR-TB-Drug-Discovery

# Install agent dependencies
pip install -r requirements_agent.txt

# Set your API key
cp .env.example .env
# Edit .env and set: ANTHROPIC_API_KEY=sk-ant-...
```

### Optional: verify install

```bash
python -c "from agent import DrugDiscoveryAgent; print('OK')"
python -m agent --help
```

---

## Usage

### Command Line

```bash
python -m agent "Find novel DprE1 inhibitors for drug-resistant tuberculosis"
python -m agent "Identify EGFR kinase inhibitors for lung cancer" --max-iterations 30
python -m agent "Discover ACE2 inhibitors for COVID-19" --output-dir my_results --quiet
```

### Python API

```python
from agent import DrugDiscoveryAgent, AgentConfig

# Basic usage
agent = DrugDiscoveryAgent(config=AgentConfig.from_env())
report = agent.run("Find InhA inhibitors for MDR-TB")

# Access results
print(report.top_compounds)    # list of dicts: {id, smiles, score, rank}
print(report.report_path)      # Path to report.md
print(report.metadata)         # dict: model, iterations, timestamp, warnings

# Custom config
cfg = AgentConfig.from_env()
cfg.model = "claude-opus-4-6"
cfg.max_iterations = 40
cfg.max_compounds_to_rank = 50
agent = DrugDiscoveryAgent(config=cfg)
report = agent.run("Find KRAS G12C inhibitors for pancreatic cancer", verbose=True)
```

---

## Tools Reference

### ChEMBL Tools

#### `search_chembl_target`

Searches ChEMBL for protein targets.

```python
result = tools.search_target(
    target_name="InhA",
    organism="Mycobacterium tuberculosis",  # optional filter
    max_results=5
)
# Returns: list of {chembl_id, pref_name, organism, target_type}
```

#### `fetch_chembl_activities`

Retrieves measured bioactivity data for a target.

```python
result = tools.fetch_activities(
    target_chembl_id="CHEMBL1849",
    activity_type="IC50",   # IC50, Ki, Kd, or EC50
    max_compounds=300
)
# Returns: list of {molecule_chembl_id, canonical_smiles, ic50_nm, pchembl_value, assay_chembl_id}
```

#### `get_chembl_compound`

Gets full molecule record.

```python
result = tools.get_compound("CHEMBL3125270")
# Returns: {canonical_smiles, molecular_formula, mw, drug_status, synonyms}
```

### PubChem Tools

#### `search_pubchem_compound`

```python
result = tools.search_compound(query="triclosan", query_type="name", max_results=5)
# query_type: "name" or "cid"
# Returns: list of {cid, canonical_smiles, molecular_formula, mw, iupac_name}
```

#### `get_pubchem_bioassays`

```python
result = tools.get_bioassays(cid=5281, max_assays=20)
# Returns: {active_count, inactive_count, assays: [...]}
```

### PDB Tools

#### `search_pdb_structure`

```python
result = tools.search_structure(
    protein_name="InhA",
    organism="Mycobacterium tuberculosis",
    max_results=5
)
# Returns: list of {pdb_id, title, resolution_angstrom, method, year}
```

#### `download_pdb_structure`

```python
result = tools.download_structure("4TZK")
# Returns: {pdb_id, resolution, chains, bound_ligands, local_path}
# File is cached at config.pdb_cache_dir / "4TZK.pdb"
```

### PubMed/EuropePMC Tools

#### `search_literature`

```python
result = tools.search_literature(
    query="InhA inhibitor tuberculosis",
    max_results=10,
    from_year=2015,
    open_access_only=False
)
# Returns: list of {pmid, title, authors, journal, year, doi, abstract}
```

#### `fetch_paper_details`

```python
result = tools.fetch_paper_details("38234567")  # PMID or DOI
# Returns: {pmid, title, authors, abstract, doi, journal, year, citations}
```

### RDKit Tools

#### `compute_admet`

Computes physicochemical and ADMET properties locally (no API call).

```python
result = tools.compute_admet(
    smiles_list=["CCc1cc(C(=O)N...)n(CC)n1", ...],
    ids=["CHEMBL3125270", ...]  # optional
)
# Returns per compound: {id, smiles, MW, cLogP, TPSA, HBD, HBA, QED,
#                        lipinski_pass, veber_pass, hERG_risk, gi_absorption, bbb_penetrant}
```

#### `train_qsar_model`

Trains a binary classifier (active/inactive) on ECFP4 fingerprints.

```python
result = tools.train_qsar_model(
    compounds=[
        {"id": "CHEMBL3125270", "smiles": "CCc1cc...", "ic50_nm": 4.0},
        ...
    ],
    activity_threshold_nm=1000.0  # IC50 < threshold → "active"
)
# Requires ≥ 30 compounds
# Returns: {best_model, roc_auc, accuracy, cv_results, feature_importance}
```

#### `rank_compounds`

Composite ranking using pIC50 + LBVS Tanimoto + QED + Lipinski.

```python
result = tools.rank_compounds(
    compounds=[{"id": ..., "smiles": ..., "ic50_nm": ...}, ...],
    reference_smiles=["Cc1ccc(Cl)cc1Cl", ...],  # known actives for LBVS
    top_n=10
)
# Returns: ranked list with composite_score, rank, and component scores
```

---

## Guardrails

The `GuardrailEnforcer` wraps every tool call and enforces data quality.

### Hard Failures (raise `GuardrailViolation`)

| Condition | Example |
|-----------|---------|
| `result.success=True` but `result.citations` is empty | Tool returned data without citing a source |
| `result.confidence` outside [0.0, 1.0] | Confidence = 1.5 or -0.1 |
| Data contains fabrication patterns | SMILES contains `"FAKE"`, confidence rationale mentions `"random"` |
| `result.success=False` but `result.error` is None/empty | Silent failure |

When a hard failure occurs, the tool returns `{"success": false, "error": "GuardrailViolation: ..."}` to Claude. Claude must handle this gracefully and try an alternative approach.

### Soft Warnings (logged, do not stop execution)

| Condition | Warning message |
|-----------|----------------|
| `result.confidence < 0.30` | "Low confidence result: ..." |
| Compound count < 50 for QSAR training | "Small dataset: QSAR results may not generalise" |

---

## MCP Server

The agent's 12 tools are also exposed as an MCP server, allowing any MCP-compatible client (Claude Desktop, Cursor, custom apps) to use them directly.

### Starting the Server

```bash
# stdio (default — for Claude Desktop)
python -m agent.mcp

# HTTP (for web clients or testing with curl)
python -m agent.mcp --http --host localhost --port 8765
```

### Testing the HTTP Server

```bash
# List all tools
curl -s http://localhost:8765/tools | python -m json.tool

# Call a tool
curl -X POST http://localhost:8765/call \
  -H "Content-Type: application/json" \
  -d '{"name": "search_chembl_target", "arguments": {"target_name": "InhA", "organism": "Mycobacterium tuberculosis"}}'
```

### Claude Desktop Configuration

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "drug-discovery": {
      "command": "python",
      "args": ["-m", "agent.mcp"],
      "cwd": "/absolute/path/to/MDR-TB-Drug-Discovery",
      "env": {
        "ANTHROPIC_API_KEY": "sk-ant-..."
      }
    }
  }
}
```

After saving, restart Claude Desktop. The 12 drug discovery tools will appear in the tools panel.

---

## Report Format

Each run saves three files to `outputs/agent_runs/<run_id>/`:

### `report.md`

Human-readable Markdown report with:
- Run metadata (prompt, model, timestamp, iterations used)
- Top compounds table (ID, SMILES, IC50, composite score)
- QSAR model performance (AUC, accuracy)
- ADMET summary (Lipinski pass rate, hERG risk flags)
- Literature citations used
- Warnings and confidence notes

### `report.json`

Machine-readable JSON with the same content, suitable for downstream processing:

```json
{
  "run_id": "20260321_143022",
  "prompt": "Find InhA inhibitors for MDR-TB",
  "model": "claude-sonnet-4-6",
  "iterations_used": 12,
  "top_compounds": [
    {
      "id": "CHEMBL3125270",
      "smiles": "CCc1cc(C(=O)N...)n(CC)n1",
      "ic50_nm": 4.0,
      "composite_score": 0.825,
      "rank": 1
    }
  ],
  "qsar_auc": 0.979,
  "warnings": [],
  "citations": [...]
}
```

### `compounds.csv`

Flat CSV of all ranked compounds with all ADMET and scoring columns, ready for analysis in R, Python, or Excel.

---

## Configuration

All settings are in `agent/config.py` and can be overridden via environment variables or the `AgentConfig` dataclass:

| Setting | Env var | Default | Description |
|---------|---------|---------|-------------|
| API key | `ANTHROPIC_API_KEY` | required | Anthropic API key |
| Model | `AGENT_MODEL` | `claude-sonnet-4-6` | Claude model ID |
| Max iterations | `AGENT_MAX_ITERATIONS` | 25 | Agentic loop cap |
| Output dir | — | `outputs/agent_runs` | Report output path |
| PDB cache | — | `data/pdb/` | Local PDB file cache |
| Max compounds | — | 500 | Max compounds fetched from ChEMBL |
| ADMET MW limit | — | 500 | Lipinski MW threshold |

---

## Extending the Agent

### Adding a New Tool

1. Create `agent/tools/mytool.py` implementing `BaseTool`:

```python
from .base import BaseTool, ToolResult, Citation

class MyTools(BaseTool):
    def my_new_tool(self, param: str) -> ToolResult:
        # ... fetch real data ...
        return ToolResult(
            tool_name="my_new_tool",
            success=True,
            data={"result": ...},
            citations=[Citation(source="MyDB", url="https://...", year=2024)],
            confidence=0.85,
            confidence_rationale="Direct database query with verified results"
        )
```

2. Export it in `agent/tools/__init__.py`
3. Add the tool dispatch in `agent/orchestrator.py` (tools dict + Claude tool schema)
4. Register it in `agent/mcp/server.py` with `@mcp.tool()`

### Swapping the Model

```python
cfg = AgentConfig.from_env()
cfg.model = "claude-opus-4-6"  # more capable, slower
agent = DrugDiscoveryAgent(config=cfg)
```

---

## Troubleshooting

### `ANTHROPIC_API_KEY` not set

```
EnvironmentError: ANTHROPIC_API_KEY is required. Set it in your environment or .env file.
```

Fix: `export ANTHROPIC_API_KEY=sk-ant-...` or add it to `.env`

### ChEMBL rate limit errors

The ChEMBL client respects a 0.4s inter-request delay. If you hit 429 errors, increase `cfg.chembl_rate_limit_delay`.

### PDB download failures

PDB files are cached at `data/pdb/`. If a download fails, delete the cached file and retry:
```bash
rm data/pdb/4TZK.pdb
```

### QSAR needs ≥ 30 compounds

The QSAR tool requires at least 30 compounds with measured IC50 values. For targets with sparse data, use `fetch_chembl_activities` with `activity_type="Ki"` or `activity_type="EC50"` to supplement.

### RDKit not installed

On some systems RDKit must be installed via conda:
```bash
conda install -c conda-forge rdkit
```
or with pip (>=2024.03.1):
```bash
pip install rdkit
```
