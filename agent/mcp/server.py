"""
MCP Server — exposes all 12 drug discovery tools via Model Context Protocol.

Usage:
    python -m agent.mcp          # stdio transport (for Claude Desktop / Cursor)
    python -m agent.mcp --http  # streamable-HTTP transport (for web MCP clients)

Claude Desktop config (claude_desktop_config.json):
    {
      "mcpServers": {
        "drug-discovery": {
          "command": "python",
          "args": ["-m", "agent.mcp"],
          "cwd": "<path to repo root>"
        }
      }
    }
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Add repo root to path so agent package is importable
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from mcp.server.fastmcp import FastMCP
from agent.tools import ChEMBLTools, PubChemTools, PDBTools, PubMedTools, RDKitTools
from agent.guardrails import GuardrailEnforcer, GuardrailViolation
from agent.config import AgentConfig

# ── Singletons ─────────────────────────────────────────────────────────────
_cfg = AgentConfig.from_env()
_chembl = ChEMBLTools()
_pubchem = PubChemTools()
_pdb = PDBTools(pdb_cache_dir=_cfg.pdb_cache_dir)
_pubmed = PubMedTools()
_rdkit = RDKitTools()
_guard = GuardrailEnforcer()

mcp = FastMCP("drug-discovery-agent")


def _wrap(result) -> str:
    """Apply guardrails and return JSON string."""
    try:
        validated = _guard.validate(result)
    except GuardrailViolation as e:
        return json.dumps({"success": False, "error": str(e)})
    return validated.to_claude_result()


# ── ChEMBL tools ──────────────────────────────────────────────────────────

@mcp.tool()
def search_chembl_target(target_name: str, organism: str = "",
                         max_results: int = 5) -> str:
    """
    Search ChEMBL for protein targets by name and optional organism.
    Returns ChEMBL target IDs, preferred names, organism, and target type.
    Use this first to identify the canonical ChEMBL ID for a protein target.

    Args:
        target_name: Protein name e.g. 'InhA', 'DprE1', 'EGFR', 'ACE2'
        organism: Organism filter e.g. 'Mycobacterium tuberculosis', 'Homo sapiens'
        max_results: Maximum results (1-20)
    """
    return _wrap(_chembl.search_target(target_name, organism, max_results))


@mcp.tool()
def fetch_chembl_activities(target_chembl_id: str,
                            activity_type: str = "IC50",
                            max_compounds: int = 300) -> str:
    """
    Retrieve measured bioactivity data from ChEMBL for a target.
    All values are from real published assays.
    Returns compound SMILES, IC50 in nM, assay ID, and pChEMBL value.

    Args:
        target_chembl_id: ChEMBL target ID e.g. 'CHEMBL1849'
        activity_type: IC50, Ki, Kd, or EC50
        max_compounds: Maximum compounds to retrieve (1-1000)
    """
    return _wrap(_chembl.fetch_activities(target_chembl_id, activity_type, max_compounds))


@mcp.tool()
def get_chembl_compound(molecule_chembl_id: str) -> str:
    """
    Get full molecule record from ChEMBL: canonical SMILES, formula, drug status, synonyms.

    Args:
        molecule_chembl_id: ChEMBL molecule ID e.g. 'CHEMBL25'
    """
    return _wrap(_chembl.get_compound(molecule_chembl_id))


# ── PubChem tools ─────────────────────────────────────────────────────────

@mcp.tool()
def search_pubchem_compound(query: str, query_type: str = "name",
                            max_results: int = 5) -> str:
    """
    Search PubChem for compounds by name or CID.
    Returns canonical SMILES, molecular formula, MW, CID, and IUPAC name.

    Args:
        query: Compound name or CID
        query_type: 'name' or 'cid'
        max_results: Maximum results (1-20)
    """
    return _wrap(_pubchem.search_compound(query, query_type, max_results))


@mcp.tool()
def get_pubchem_bioassays(cid: int, max_assays: int = 20) -> str:
    """
    Retrieve PubChem BioAssay activity summary for a compound CID.
    Returns active/inactive assay counts and activity outcomes.

    Args:
        cid: PubChem compound CID (numeric identifier)
        max_assays: Maximum assays to return (1-100)
    """
    return _wrap(_pubchem.get_bioassays(cid, max_assays))


# ── PDB tools ─────────────────────────────────────────────────────────────

@mcp.tool()
def search_pdb_structure(protein_name: str, organism: str = "",
                         max_results: int = 5) -> str:
    """
    Search RCSB PDB for protein structures. Returns PDB IDs, resolution, method.
    Prefer structures with resolution < 2.5 Å for drug discovery.

    Args:
        protein_name: Protein name e.g. 'InhA', 'DprE1', 'EGFR'
        organism: Scientific name filter e.g. 'Mycobacterium tuberculosis'
        max_results: Maximum structures (1-20)
    """
    return _wrap(_pdb.search_structure(protein_name, organism, max_results))


@mcp.tool()
def download_pdb_structure(pdb_id: str) -> str:
    """
    Download a PDB structure from RCSB and extract metadata.
    Returns resolution, chains, bound ligands. File cached locally.

    Args:
        pdb_id: 4-character PDB entry ID e.g. '4TZK', '1IVA'
    """
    return _wrap(_pdb.download_structure(pdb_id))


# ── PubMed/EuropePMC tools ────────────────────────────────────────────────

@mcp.tool()
def search_literature(query: str, max_results: int = 10,
                      from_year: int = 2010,
                      open_access_only: bool = False) -> str:
    """
    Search EuropePMC/PubMed for peer-reviewed publications.
    Returns titles, authors, journal, year, DOIs, and abstracts.

    Args:
        query: Search terms e.g. 'InhA inhibitor tuberculosis'
        max_results: Maximum papers (1-50)
        from_year: Filter to papers from this year onward
        open_access_only: Return only open-access papers
    """
    return _wrap(_pubmed.search_literature(query, max_results, from_year, open_access_only))


@mcp.tool()
def fetch_paper_details(pmid: str) -> str:
    """
    Fetch full abstract and metadata for a paper by PubMed ID or DOI.

    Args:
        pmid: PubMed ID (numeric) or DOI string
    """
    return _wrap(_pubmed.fetch_paper_details(pmid))


# ── RDKit tools ───────────────────────────────────────────────────────────

@mcp.tool()
def compute_admet(smiles_list: list[str], ids: list[str] | None = None) -> str:
    """
    Compute ADMET properties for SMILES strings using RDKit (local, no API call).
    Returns MW, cLogP, TPSA, HBD, HBA, QED, Lipinski/Veber compliance, hERG flag.

    Args:
        smiles_list: List of SMILES strings (max 500)
        ids: Optional parallel list of compound IDs
    """
    return _wrap(_rdkit.compute_admet(smiles_list, ids))


@mcp.tool()
def train_qsar_model(compounds: list[dict],
                     activity_threshold_nm: float = 1000.0) -> str:
    """
    Train a QSAR binary classifier on compounds with measured IC50 values.
    Evaluates Random Forest, Gradient Boosting, Logistic Regression with 5-fold CV.
    Requires ≥30 compounds with SMILES and IC50.

    Args:
        compounds: List of {id, smiles, ic50_nm} dicts
        activity_threshold_nm: IC50 below this (nM) is 'active' (default 1000)
    """
    return _wrap(_rdkit.train_qsar_model(compounds, activity_threshold_nm))


@mcp.tool()
def rank_compounds(compounds: list[dict],
                   reference_smiles: list[str] | None = None,
                   top_n: int = 10) -> str:
    """
    Rank compounds by composite score: pIC50 (40%) + LBVS Tanimoto (30%) + QED (20%) + Lipinski (10%).
    All IC50 values must be real measured values. Returns ranked table.

    Args:
        compounds: List of {id, smiles, ic50_nm} dicts
        reference_smiles: Known active reference SMILES for LBVS scoring
        top_n: Number of top compounds to return (1-50)
    """
    return _wrap(_rdkit.rank_compounds(compounds, reference_smiles, top_n))


# ── Entry point ───────────────────────────────────────────────────────────

def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Drug Discovery MCP Server")
    parser.add_argument("--http", action="store_true",
                        help="Use streamable-HTTP transport instead of stdio")
    parser.add_argument("--host", default="localhost",
                        help="HTTP host (default: localhost)")
    parser.add_argument("--port", type=int, default=8765,
                        help="HTTP port (default: 8765)")
    args = parser.parse_args()

    if args.http:
        print(f"[MCP] Starting HTTP server at http://{args.host}:{args.port}", file=sys.stderr)
        mcp.run(transport="streamable-http", host=args.host, port=args.port)
    else:
        mcp.run()  # stdio — default for Claude Desktop


if __name__ == "__main__":
    main()
