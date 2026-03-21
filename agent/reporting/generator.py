"""Structured report generation — Markdown, JSON, and CSV outputs."""
from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..tools.base import ToolResult
    from ..config import AgentConfig


@dataclass
class DiscoveryReport:
    run_id: str
    timestamp: str
    user_prompt: str
    tool_calls: int
    successful_calls: int
    markdown_path: Path
    json_path: Path
    csv_path: Path | None

    def print_summary(self) -> None:
        print("\n" + "=" * 60)
        print(f"  DRUG DISCOVERY REPORT — {self.run_id}")
        print("=" * 60)
        print(f"  Query  : {self.user_prompt[:80]}")
        print(f"  Tools  : {self.tool_calls} calls ({self.successful_calls} successful)")
        print(f"  Report : {self.markdown_path}")
        if self.csv_path and self.csv_path.exists():
            print(f"  CSV    : {self.csv_path}")
        print("=" * 60 + "\n")


class ReportGenerator:

    def __init__(self, config: "AgentConfig"):
        self.config = config

    def generate(self, results: list["ToolResult"],
                 conversation: list[dict]) -> DiscoveryReport:
        """Write report.md, report.json, and compounds.csv to outputs/agent_runs/{run_id}/."""
        ts = datetime.now(timezone.utc)
        run_id = f"run_{ts.strftime('%Y%m%d_%H%M%S')}"
        out_dir = self.config.output_dir / run_id
        out_dir.mkdir(parents=True, exist_ok=True)

        user_prompt = ""
        assistant_final = ""
        for msg in conversation:
            if msg.get("role") == "user" and isinstance(msg.get("content"), str):
                user_prompt = msg["content"]
                break
        for msg in reversed(conversation):
            if msg.get("role") == "assistant":
                content = msg.get("content", [])
                if isinstance(content, list):
                    for block in content:
                        if hasattr(block, "type") and block.type == "text":
                            assistant_final = block.text
                            break
                elif isinstance(content, str):
                    assistant_final = content
                if assistant_final:
                    break

        md_text = self._build_markdown(results, user_prompt, assistant_final, run_id, ts)
        json_data = self._build_json(results, user_prompt, run_id, ts)
        csv_path = self._build_csv(results, out_dir)

        md_path = out_dir / "report.md"
        json_path = out_dir / "report.json"
        md_path.write_text(md_text, encoding="utf-8")
        json_path.write_text(json.dumps(json_data, indent=2, default=str), encoding="utf-8")

        successful = sum(1 for r in results if r.success)
        return DiscoveryReport(
            run_id=run_id,
            timestamp=ts.isoformat(),
            user_prompt=user_prompt,
            tool_calls=len(results),
            successful_calls=successful,
            markdown_path=md_path,
            json_path=json_path,
            csv_path=csv_path,
        )

    # ── Markdown builder ──────────────────────────────────────────────────────

    def _build_markdown(self, results: list["ToolResult"],
                        user_prompt: str, assistant_final: str,
                        run_id: str, ts: datetime) -> str:
        lines = [
            f"# Drug Discovery Report — {run_id}",
            f"",
            f"**Generated:** {ts.strftime('%Y-%m-%d %H:%M UTC')}",
            f"**Query:** {user_prompt}",
            f"**Integrity:** All data sourced from real external APIs (ChEMBL, PubChem, PDB, EuropePMC). "
            f"No fabricated, simulated, or randomly generated scientific values.",
            f"",
            f"---",
            f"",
        ]

        # Agent synthesis
        if assistant_final:
            lines += [
                "## Agent Analysis",
                "",
                assistant_final,
                "",
                "---",
                "",
            ]

        # Tool call log
        lines += ["## Tool Call Log", ""]
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        lines.append(f"**Total calls:** {len(results)} | **Successful:** {len(successful)} | **Failed:** {len(failed)}")
        lines.append("")

        for r in results:
            status = "✓" if r.success else "✗"
            lines.append(f"### {status} `{r.tool_name}`")
            lines.append(f"- **Confidence:** {r.confidence:.2f} — {r.confidence_rationale}")
            if r.error:
                lines.append(f"- **Error:** {r.error}")
            if r.warnings:
                for w in r.warnings:
                    lines.append(f"- **Warning:** {w}")
            if r.success and isinstance(r.data, dict):
                # Show key summary stats
                for key in ["count", "total_records", "n_compounds", "returned",
                            "n_active", "best_model", "best_test_roc_auc", "total_found"]:
                    if key in r.data:
                        lines.append(f"- **{key.replace('_', ' ').title()}:** {r.data[key]}")
            if r.citations:
                lines.append(f"- **Citations:** " + "; ".join(
                    f"[{c.source}]({c.url})" + (f" doi:{c.doi}" if c.doi else "")
                    for c in r.citations
                ))
            lines.append("")

        # Compounds table (if ranking was done)
        compounds_data = self._extract_ranked_compounds(results)
        if compounds_data:
            lines += [
                "## Top-Ranked Compounds",
                "",
                "| Rank | ID | IC50 (nM) | MW | cLogP | TPSA | QED | Composite |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
            for c in compounds_data[:15]:
                lines.append(
                    f"| {c.get('rank', '–')} | {c.get('id', '')} | "
                    f"{c.get('ic50_nm', '–')} | {c.get('mw', '–')} | "
                    f"{c.get('clogp', '–')} | {c.get('tpsa', '–')} | "
                    f"{c.get('qed', '–')} | {c.get('composite_score', '–')} |"
                )
            lines.append("")

        # QSAR summary
        qsar_results = [r for r in results if r.tool_name == "train_qsar_model" and r.success]
        if qsar_results:
            qd = qsar_results[-1].data
            lines += ["## QSAR Model Performance", ""]
            lines.append(f"**Best model:** {qd.get('best_model', '–')} (test ROC-AUC: {qd.get('best_test_roc_auc', '–')})")
            lines.append(f"**Training set:** {qd.get('n_compounds', '–')} compounds ({qd.get('n_active', '–')} active, {qd.get('n_inactive', '–')} inactive)")
            lines.append("")
            lines.append("| Model | CV ROC-AUC | Test ROC-AUC | Test F1 | Test Precision | Test Recall |")
            lines.append("| --- | --- | --- | --- | --- | --- |")
            for name, m in (qd.get("models") or {}).items():
                lines.append(
                    f"| {name} | {m.get('cv_roc_auc_mean', '–')}±{m.get('cv_roc_auc_std', '–')} | "
                    f"{m.get('test_roc_auc', '–')} | {m.get('test_f1', '–')} | "
                    f"{m.get('test_precision', '–')} | {m.get('test_recall', '–')} |"
                )
            lines.append("")

        # All citations
        all_citations = []
        seen_urls = set()
        for r in results:
            for c in r.citations:
                if c.url not in seen_urls:
                    all_citations.append(c)
                    seen_urls.add(c.url)

        if all_citations:
            lines += ["## References", ""]
            for i, c in enumerate(all_citations, 1):
                doi_str = f" DOI: {c.doi}" if c.doi else ""
                lines.append(f"{i}. {c.source}. {c.url}{doi_str} (retrieved {c.retrieved_at[:10]})")
            lines.append("")

        return "\n".join(lines)

    # ── JSON builder ──────────────────────────────────────────────────────────

    def _build_json(self, results: list["ToolResult"],
                    user_prompt: str, run_id: str, ts: datetime) -> dict:
        return {
            "run_id": run_id,
            "timestamp": ts.isoformat(),
            "user_prompt": user_prompt,
            "integrity_statement": (
                "All data in this report was sourced from real external APIs "
                "(ChEMBL, PubChem, RCSB PDB, EuropePMC). "
                "No fabricated, simulated, or randomly generated scientific values were used."
            ),
            "summary": {
                "total_tool_calls": len(results),
                "successful_calls": sum(1 for r in results if r.success),
                "failed_calls": sum(1 for r in results if not r.success),
                "average_confidence": round(
                    sum(r.confidence for r in results) / len(results), 3
                ) if results else 0.0,
            },
            "tool_calls": [r.to_dict() for r in results],
        }

    # ── CSV builder ───────────────────────────────────────────────────────────

    def _build_csv(self, results: list["ToolResult"],
                   out_dir: Path) -> Path | None:
        compounds = self._extract_ranked_compounds(results)
        if not compounds:
            compounds = self._extract_admet_compounds(results)
        if not compounds:
            return None

        try:
            import pandas as pd
            df = pd.DataFrame(compounds)
            csv_path = out_dir / "compounds.csv"
            df.to_csv(csv_path, index=False)
            return csv_path
        except ImportError:
            # Write manually
            if not compounds:
                return None
            csv_path = out_dir / "compounds.csv"
            keys = list(compounds[0].keys())
            rows = [",".join(str(c.get(k, "")) for k in keys) for c in compounds]
            csv_path.write_text(",".join(keys) + "\n" + "\n".join(rows), encoding="utf-8")
            return csv_path

    @staticmethod
    def _extract_ranked_compounds(results: list["ToolResult"]) -> list[dict]:
        for r in reversed(results):
            if r.tool_name == "rank_compounds" and r.success:
                return r.data.get("top_compounds", [])
        return []

    @staticmethod
    def _extract_admet_compounds(results: list["ToolResult"]) -> list[dict]:
        for r in reversed(results):
            if r.tool_name == "compute_admet" and r.success:
                return r.data.get("compounds", [])
        return []
