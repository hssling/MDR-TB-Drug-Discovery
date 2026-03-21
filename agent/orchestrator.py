"""
Core agentic loop — Claude tool_use orchestration with guardrails.

Flow:
  1. User prompt → Claude with all 15 tool definitions
  2. Claude calls tools (tool_use blocks)
  3. Each call dispatched → real API/computation → GuardrailEnforcer
  4. Results fed back as tool_result messages
  5. Loop until stop_reason == "end_turn"
  6. ReportGenerator writes Markdown/JSON/CSV
"""
from __future__ import annotations

import sys
from typing import Callable
from .config import AgentConfig
from .tools import ChEMBLTools, PubChemTools, PDBTools, PubMedTools, RDKitTools, ToolResult
from .guardrails import GuardrailEnforcer, GuardrailViolation
from .reporting import ReportGenerator, DiscoveryReport

_SYSTEM_PROMPT = """You are a rigorous, evidence-based computational drug discovery scientist.

Your mission: given a disease or protein target from the user, conduct a complete, structured drug discovery analysis using ONLY the provided tools.

MANDATORY RULES — violations are not tolerated:
1. ALL scientific claims must be supported by tool results. Do NOT invent numbers, compounds, mechanisms, or statistics.
2. If a tool fails, report the failure as a limitation. Never substitute invented data.
3. Every compound IC50 you report must come from a tool result, not from memory.
4. Always include citations from tool results in your final report.
5. Report confidence scores and their rationale honestly.
6. When QSAR is possible (≥30 compounds with IC50), always train and report it.
7. Your final answer must be a complete structured report with all sections below.

REQUIRED REPORT STRUCTURE:
## 1. Disease & Target Context
## 2. Target Identification (from ChEMBL/PubMed)
## 3. Structural Context (PDB entry, resolution, key binding residues)
## 4. Compound Dataset (ChEMBL bioactivity summary)
## 5. QSAR Model Performance (if dataset ≥30 compounds)
## 6. ADMET Analysis (Lipinski/Veber/hERG summary)
## 7. Top-Ranked Compounds (composite score table, all real IC50 values)
## 8. Lead Compound Profile (top hit: SMILES, ADMET, structural rationale)
## 9. Literature Support (key papers from EuropePMC/PubMed)
## 10. Limitations & Next Steps (be specific about what was NOT done)

WORKFLOW GUIDANCE:
Step 1: search_literature for disease context and known targets
Step 2: search_chembl_target to identify the canonical ChEMBL target ID
Step 3: search_pdb_structure + download_pdb_structure for structural context
Step 4: fetch_chembl_activities to get real IC50 data
Step 5: compute_admet on the active compounds
Step 6: train_qsar_model if ≥30 compounds available
Step 7: search_pubchem_compound for any reference compounds you want to look up
Step 8: rank_compounds to produce the final prioritized list
Step 9: search_literature for specific compound series and mechanism papers
Step 10: synthesize all results into the required report structure"""


class DrugDiscoveryOrchestrator:
    """Main agent class. Call run(prompt) to execute a full discovery workflow."""

    def __init__(self, config: AgentConfig | None = None):
        self.config = config or AgentConfig.from_env()
        self.config.validate()

        # Tool instances
        self._chembl = ChEMBLTools()
        self._pubchem = PubChemTools()
        self._pdb = PDBTools(pdb_cache_dir=self.config.pdb_cache_dir)
        self._pubmed = PubMedTools()
        self._rdkit = RDKitTools()
        self._enforcer = GuardrailEnforcer()

        # Aggregate tool definitions for Claude
        self._tool_defs: list[dict] = (
            self._chembl.get_tool_definitions() +
            self._pubchem.get_tool_definitions() +
            self._pdb.get_tool_definitions() +
            self._pubmed.get_tool_definitions() +
            self._rdkit.get_tool_definitions()
        )

        # Dispatch table
        self._dispatch: dict[str, Callable] = {
            "search_chembl_target":         lambda **kw: self._chembl.search_target(**kw),
            "fetch_chembl_activities":       lambda **kw: self._chembl.fetch_activities(**kw),
            "get_chembl_compound":           lambda **kw: self._chembl.get_compound(**kw),
            "search_pubchem_compound":       lambda **kw: self._pubchem.search_compound(**kw),
            "get_pubchem_bioassays":         lambda **kw: self._pubchem.get_bioassays(**kw),
            "search_pdb_structure":          lambda **kw: self._pdb.search_structure(**kw),
            "download_pdb_structure":        lambda **kw: self._pdb.download_structure(**kw),
            "search_literature":             lambda **kw: self._pubmed.search_literature(**kw),
            "fetch_paper_details":           lambda **kw: self._pubmed.fetch_paper_details(**kw),
            "compute_admet":                 lambda **kw: self._rdkit.compute_admet(**kw),
            "train_qsar_model":              lambda **kw: self._rdkit.train_qsar_model(**kw),
            "rank_compounds":                lambda **kw: self._rdkit.rank_compounds(**kw),
        }

        self._all_results: list[ToolResult] = []
        self._messages: list[dict] = []

    # ── Public API ─────────────────────────────────────────────────────────────

    def run(self, user_prompt: str, verbose: bool = True) -> DiscoveryReport:
        """Execute a complete drug discovery workflow for the given prompt."""
        import anthropic

        client = anthropic.Anthropic(api_key=self.config.anthropic_api_key)
        self._all_results = []
        self._messages = [{"role": "user", "content": user_prompt}]

        if verbose:
            print(f"\n[Agent] Starting: {user_prompt[:80]}...")
            print(f"[Agent] Model: {self.config.model} | Max iterations: {self.config.max_iterations}")

        for iteration in range(self.config.max_iterations):
            response = client.messages.create(
                model=self.config.model,
                max_tokens=self.config.max_tokens,
                system=_SYSTEM_PROMPT,
                tools=self._tool_defs,
                messages=self._messages,
            )

            # Add assistant turn to history
            self._messages.append({
                "role": "assistant",
                "content": response.content
            })

            if verbose:
                self._log_response(iteration, response)

            if response.stop_reason == "end_turn":
                break

            if response.stop_reason == "tool_use":
                tool_results_content = []
                for block in response.content:
                    if block.type == "tool_use":
                        result = self._execute_tool(block.name, block.input, verbose=verbose)
                        self._all_results.append(result)
                        tool_results_content.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": result.to_claude_result(),
                        })

                self._messages.append({
                    "role": "user",
                    "content": tool_results_content,
                })
            else:
                if verbose:
                    print(f"[Agent] Unexpected stop_reason: {response.stop_reason}")
                break

        if verbose:
            print(f"\n[Agent] Generating report ({len(self._all_results)} tool calls)...")

        report = ReportGenerator(self.config).generate(
            results=self._all_results,
            conversation=self._messages,
        )
        return report

    # ── Internal helpers ───────────────────────────────────────────────────────

    def _execute_tool(self, tool_name: str, tool_input: dict,
                      verbose: bool = True) -> ToolResult:
        if verbose:
            print(f"  → {tool_name}({', '.join(f'{k}={repr(v)[:40]}' for k, v in tool_input.items())})")

        if tool_name not in self._dispatch:
            return ToolResult(
                tool_name=tool_name,
                success=False,
                data={},
                citations=[],
                confidence=0.0,
                confidence_rationale="Unknown tool",
                error=f"Tool '{tool_name}' is not registered in the dispatch table.",
            )

        try:
            result = self._dispatch[tool_name](**tool_input)
            validated = self._enforcer.validate(result)
            if verbose and not validated.success:
                print(f"    ✗ FAILED: {validated.error}")
            elif verbose:
                print(f"    ✓ confidence={validated.confidence:.2f}")
            return validated

        except GuardrailViolation as e:
            if verbose:
                print(f"    ✗ GUARDRAIL VIOLATION: {e}")
            return ToolResult(
                tool_name=tool_name,
                success=False,
                data={"guardrail_violation": str(e)},
                citations=[],
                confidence=0.0,
                confidence_rationale="Blocked by guardrail enforcer",
                error=f"GuardrailViolation: {e}",
            )

        except Exception as exc:
            if verbose:
                print(f"    ✗ ERROR: {exc}")
            return ToolResult(
                tool_name=tool_name,
                success=False,
                data={},
                citations=[],
                confidence=0.0,
                confidence_rationale="Tool execution raised an unhandled exception",
                error=str(exc),
            )

    @staticmethod
    def _log_response(iteration: int, response) -> None:
        tool_names = [b.name for b in response.content if hasattr(b, "name")]
        text_blocks = [b for b in response.content if hasattr(b, "text")]
        if tool_names:
            print(f"\n[Iter {iteration+1}] Claude calling tools: {tool_names}")
        elif text_blocks:
            preview = text_blocks[0].text[:120].replace("\n", " ")
            print(f"\n[Iter {iteration+1}] Claude response (stop={response.stop_reason}): {preview}...")
        else:
            print(f"\n[Iter {iteration+1}] stop_reason={response.stop_reason}")
