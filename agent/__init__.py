"""
Drug Discovery Agent
====================
Disease-agnostic, Claude-orchestrated, real-data-only computational drug discovery.

Usage:
    python -m agent "Find novel DprE1 inhibitors for drug-resistant tuberculosis"
    python -m agent.mcp  # start MCP server for Claude Desktop / Cursor integration

API:
    from agent import DrugDiscoveryAgent
    report = DrugDiscoveryAgent().run("Find EGFR inhibitors for lung cancer")
"""
from .orchestrator import DrugDiscoveryOrchestrator as DrugDiscoveryAgent
from .config import AgentConfig
from .reporting.generator import DiscoveryReport

__version__ = "1.0.0"
__all__ = ["DrugDiscoveryAgent", "AgentConfig", "DiscoveryReport"]
