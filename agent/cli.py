"""CLI entry point for the Drug Discovery Agent."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="drug-discovery-agent",
        description="AI-powered, evidence-based computational drug discovery for any disease.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m agent "Find novel DprE1 inhibitors for drug-resistant tuberculosis"
  python -m agent "Identify EGFR kinase inhibitors for lung cancer"
  python -m agent "Discover ACE2 inhibitors for COVID-19"
  python -m agent "Find PfDHFR inhibitors for malaria"

MCP server:
  python -m agent.mcp                   # stdio (Claude Desktop / Cursor)
  python -m agent.mcp --http            # HTTP on localhost:8765
  python -m agent.mcp --http --port 9000
        """
    )
    parser.add_argument(
        "prompt",
        nargs="?",
        help="Drug discovery task e.g. 'Find InhA inhibitors for MDR-TB'"
    )
    parser.add_argument(
        "--model",
        default=os.environ.get("AGENT_MODEL", "claude-sonnet-4-6"),
        help="Claude model ID (default: claude-sonnet-4-6)"
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=25,
        help="Maximum agentic loop iterations (default: 25)"
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/agent_runs",
        help="Directory for report output (default: outputs/agent_runs)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress per-tool progress output"
    )
    args = parser.parse_args()

    if not args.prompt:
        parser.print_help()
        sys.exit(0)

    # Late import so --help works without ANTHROPIC_API_KEY
    from .config import AgentConfig
    from .orchestrator import DrugDiscoveryOrchestrator

    cfg = AgentConfig.from_env()
    cfg.model = args.model
    cfg.max_iterations = args.max_iterations
    cfg.output_dir = Path(args.output_dir)

    try:
        cfg.validate()
    except EnvironmentError as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)

    agent = DrugDiscoveryOrchestrator(config=cfg)
    report = agent.run(args.prompt, verbose=not args.quiet)
    report.print_summary()
