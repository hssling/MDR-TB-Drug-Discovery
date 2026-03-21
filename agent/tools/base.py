"""
Base contracts for all agent tools.

Every tool MUST return a ToolResult with:
  - citations: non-empty list (unless success=False)
  - confidence: float in [0.0, 1.0]
  - confidence_rationale: non-empty string
  - data: real data (no np.random, no mock values)
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from abc import ABC, abstractmethod


@dataclass
class Citation:
    source: str          # Human-readable name e.g. "ChEMBL database (EMBL-EBI)"
    url: str             # Canonical URL — must be non-empty
    retrieved_at: str    # ISO-8601 timestamp
    doi: str = ""        # Optional DOI

    @classmethod
    def now(cls, source: str, url: str, doi: str = "") -> "Citation":
        return cls(
            source=source,
            url=url,
            retrieved_at=datetime.now(timezone.utc).isoformat(),
            doi=doi,
        )

    def to_dict(self) -> dict:
        return {
            "source": self.source,
            "url": self.url,
            "retrieved_at": self.retrieved_at,
            "doi": self.doi,
        }


@dataclass
class ToolResult:
    tool_name: str
    success: bool
    data: dict | list
    citations: list[Citation]
    confidence: float               # [0.0, 1.0]
    confidence_rationale: str
    error: str | None = None
    warnings: list[str] = field(default_factory=list)

    def to_claude_result(self) -> str:
        """Serialize for use as Claude tool_result content (JSON string)."""
        payload = {
            "tool_name": self.tool_name,
            "success": self.success,
            "confidence": round(self.confidence, 3),
            "confidence_rationale": self.confidence_rationale,
            "citations": [c.to_dict() for c in self.citations],
            "warnings": self.warnings,
        }
        if self.error:
            payload["error"] = self.error
        else:
            payload["data"] = self.data
        return json.dumps(payload, default=str, ensure_ascii=False)

    def to_dict(self) -> dict:
        return {
            "tool_name": self.tool_name,
            "success": self.success,
            "confidence": self.confidence,
            "confidence_rationale": self.confidence_rationale,
            "citations": [c.to_dict() for c in self.citations],
            "warnings": self.warnings,
            "error": self.error,
            "data": self.data,
        }


class BaseTool(ABC):
    """Abstract base class. Subclasses implement individual tool methods."""

    @staticmethod
    def _cite(source: str, url: str, doi: str = "") -> Citation:
        return Citation.now(source=source, url=url, doi=doi)

    @staticmethod
    def _ok(tool_name: str, data: dict | list, citations: list[Citation],
            confidence: float, rationale: str,
            warnings: list[str] | None = None) -> ToolResult:
        return ToolResult(
            tool_name=tool_name,
            success=True,
            data=data,
            citations=citations,
            confidence=confidence,
            confidence_rationale=rationale,
            warnings=warnings or [],
        )

    @staticmethod
    def _err(tool_name: str, error: str) -> ToolResult:
        return ToolResult(
            tool_name=tool_name,
            success=False,
            data={},
            citations=[],
            confidence=0.0,
            confidence_rationale="Tool execution failed; no data to assess.",
            error=error,
        )

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        """Return Anthropic-format tool definitions (input_schema style)."""
        raise NotImplementedError
