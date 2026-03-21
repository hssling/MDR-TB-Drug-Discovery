"""Agent configuration — reads from environment variables, no hardcoded secrets."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class AgentConfig:
    # ── Anthropic API ──────────────────────────────────────────────────────────
    anthropic_api_key: str = field(default_factory=lambda: os.environ.get("ANTHROPIC_API_KEY", ""))
    model: str = "claude-sonnet-4-6"
    max_tokens: int = 8096
    max_iterations: int = 25          # hard cap on agentic loop iterations

    # ── External API rate limits (seconds between requests) ───────────────────
    chembl_delay: float = 0.4
    pubchem_delay: float = 0.3
    pdb_delay: float = 0.5
    europepmc_delay: float = 0.4

    # ── QSAR / ML settings ────────────────────────────────────────────────────
    qsar_test_size: float = 0.20
    qsar_cv_folds: int = 5
    qsar_random_state: int = 42
    activity_threshold_nm: float = 1000.0   # IC50 < this → "active"
    min_dataset_size: int = 30              # refuse QSAR below this

    # ── ADMET hard filters ────────────────────────────────────────────────────
    mw_min: float = 150.0
    mw_max: float = 600.0
    clogp_max: float = 5.5
    tpsa_max: float = 140.0
    hbd_max: int = 5
    hba_max: int = 10
    rotbonds_max: int = 10

    # ── Composite ranking weights ─────────────────────────────────────────────
    weight_pic50: float = 0.40
    weight_lbvs: float = 0.30
    weight_qed: float = 0.20
    weight_lipinski: float = 0.10

    # ── Paths ─────────────────────────────────────────────────────────────────
    repo_root: Path = field(
        default_factory=lambda: Path(__file__).resolve().parent.parent
    )
    output_dir: Path = field(
        default_factory=lambda: Path(__file__).resolve().parent.parent / "outputs" / "agent_runs"
    )
    pdb_cache_dir: Path = field(
        default_factory=lambda: Path(__file__).resolve().parent.parent / "data" / "pdb"
    )

    # ── MCP server ────────────────────────────────────────────────────────────
    mcp_host: str = "localhost"
    mcp_port: int = 8765

    # ── Guardrails ────────────────────────────────────────────────────────────
    require_citations: bool = True
    min_confidence_warn: float = 0.30   # warn (not fail) below this
    block_fabrication: bool = True

    @classmethod
    def from_env(cls) -> "AgentConfig":
        cfg = cls()
        cfg.anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
        cfg.model = os.environ.get("AGENT_MODEL", cfg.model)
        cfg.output_dir = Path(os.environ.get("AGENT_OUTPUT_DIR", str(cfg.output_dir)))
        return cfg

    def validate(self) -> None:
        if not self.anthropic_api_key:
            raise EnvironmentError(
                "ANTHROPIC_API_KEY environment variable is not set. "
                "Export it before running: export ANTHROPIC_API_KEY=sk-ant-..."
            )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pdb_cache_dir.mkdir(parents=True, exist_ok=True)
