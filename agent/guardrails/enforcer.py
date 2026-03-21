"""
Guardrail enforcer — intercepts every ToolResult before it reaches Claude or MCP clients.

Blocked patterns (raises GuardrailViolation):
  - Missing citations on a successful result
  - Confidence outside [0.0, 1.0]
  - Known fabrication markers in data (np.random, "mock_", "synthetic_data")
  - success=False with empty error message

Soft warnings (adds to result.warnings):
  - Confidence below 0.30
  - Citations with placeholder URLs
  - Dataset too small for QSAR
"""
from __future__ import annotations

import json
import re
from ..tools.base import ToolResult


class GuardrailViolation(Exception):
    """Raised when a tool result violates the no-fabrication contract."""
    def __init__(self, tool_name: str, violation: str):
        self.tool_name = tool_name
        self.violation = violation
        super().__init__(f"[GUARDRAIL:{tool_name}] {violation}")


class GuardrailEnforcer:
    """
    Call enforcer.validate(result) after every tool execution.
    Returns the result unchanged (possibly with added warnings),
    or raises GuardrailViolation for critical integrity failures.
    """

    # Patterns that indicate fabricated/synthetic data leaking through
    _FABRICATION_PATTERNS = [
        (r"\bnp\.random\b", "numpy random data generation detected"),
        (r"\brandom\.uniform\b", "random.uniform in data — fabricated values suspected"),
        (r'"mock_', "mock_ prefix in data keys — test/mock data in production result"),
        (r'"synthetic_data"', "synthetic_data key in result"),
        (r'"generated_randomly"', "generated_randomly flag in result"),
        (r'"placeholder"', "placeholder value in result"),
        (r"FAKE_IC50", "FAKE_IC50 marker in result"),
        (r"FABRICATED", "FABRICATED marker in result"),
    ]

    _PLACEHOLDER_URLS = {
        "", "http://example.com", "https://example.org",
        "http://placeholder.invalid", "https://todo.fix"
    }

    def validate(self, result: ToolResult) -> ToolResult:
        """
        Validate a ToolResult. Returns validated result or raises GuardrailViolation.
        May append warnings to result.warnings for soft issues.
        """
        self._check_error_state(result)
        self._check_confidence(result)

        if result.success:
            self._check_citations(result)
            self._check_no_fabrication(result)

        self._soft_checks(result)
        return result

    # ── Hard checks (raise on failure) ────────────────────────────────────────

    def _check_error_state(self, result: ToolResult) -> None:
        if not result.success and not result.error:
            raise GuardrailViolation(
                result.tool_name,
                "success=False but error field is empty — tool must provide a descriptive error message."
            )

    def _check_confidence(self, result: ToolResult) -> None:
        c = result.confidence
        if not (0.0 <= c <= 1.0):
            raise GuardrailViolation(
                result.tool_name,
                f"confidence={c} is outside [0.0, 1.0] — invalid value."
            )
        if not result.confidence_rationale:
            raise GuardrailViolation(
                result.tool_name,
                "confidence_rationale is empty — every tool must explain its confidence score."
            )

    def _check_citations(self, result: ToolResult) -> None:
        if not result.citations:
            raise GuardrailViolation(
                result.tool_name,
                "citations list is empty on a successful result. "
                "Every successful tool call must cite at least one data source."
            )
        for cite in result.citations:
            if not cite.url or cite.url in self._PLACEHOLDER_URLS:
                raise GuardrailViolation(
                    result.tool_name,
                    f"Citation has invalid/placeholder URL: '{cite.url}'. "
                    "All citations must have real, accessible URLs."
                )

    def _check_no_fabrication(self, result: ToolResult) -> None:
        try:
            data_str = json.dumps(result.data, default=str)
        except Exception:
            data_str = str(result.data)

        for pattern, message in self._FABRICATION_PATTERNS:
            if re.search(pattern, data_str):
                raise GuardrailViolation(
                    result.tool_name,
                    f"Fabrication pattern detected in result data: {message}. "
                    "Tool results must contain only real, sourced data."
                )

    # ── Soft checks (add warnings, do not raise) ──────────────────────────────

    def _soft_checks(self, result: ToolResult) -> None:
        if result.confidence < 0.30 and result.success:
            result.warnings.append(
                f"Low confidence ({result.confidence:.2f}): {result.confidence_rationale}"
            )

        if result.success and isinstance(result.data, dict):
            n = result.data.get("total_records") or result.data.get("n_compounds") or 0
            if 0 < int(n) < 30 and result.tool_name in {"fetch_chembl_activities", "train_qsar_model"}:
                result.warnings.append(
                    f"Small dataset ({n} compounds) — QSAR/ranking results should be interpreted with caution."
                )

        if result.success and result.citations:
            for cite in result.citations:
                if not cite.source or cite.source.strip() in {"", "Unknown"}:
                    result.warnings.append(
                        f"Citation has empty source name — please verify the data origin."
                    )
