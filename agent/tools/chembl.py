"""ChEMBL tools — target search, bioactivity retrieval, compound lookup."""
from __future__ import annotations

import time
from .base import BaseTool, ToolResult, Citation


class ChEMBLTools(BaseTool):

    _CHEMBL_URL = "https://www.ebi.ac.uk/chembl"
    _last_call: float = 0.0
    _delay: float = 0.4

    # ── Tool definitions for Claude ───────────────────────────────────────────

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        return [
            {
                "name": "search_chembl_target",
                "description": (
                    "Search ChEMBL for protein targets by name and optional organism. "
                    "Returns ChEMBL target IDs, preferred names, organism, and target type. "
                    "Use this first to identify the canonical ChEMBL ID for a protein target."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "target_name": {
                            "type": "string",
                            "description": "Protein name e.g. 'InhA', 'DprE1', 'EGFR', 'ACE2'"
                        },
                        "organism": {
                            "type": "string",
                            "description": "Organism filter e.g. 'Mycobacterium tuberculosis', 'Homo sapiens'. Leave empty to search all organisms."
                        },
                        "max_results": {
                            "type": "integer",
                            "description": "Maximum results to return (1-20)",
                            "default": 5,
                            "minimum": 1,
                            "maximum": 20
                        }
                    },
                    "required": ["target_name"]
                }
            },
            {
                "name": "fetch_chembl_activities",
                "description": (
                    "Retrieve measured bioactivity data (IC50, Ki, Kd, or EC50) from ChEMBL "
                    "for a given target ChEMBL ID. All values are from real published assays. "
                    "Returns compound SMILES, IC50 in nM, assay ID, and pChEMBL value."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "target_chembl_id": {
                            "type": "string",
                            "description": "ChEMBL target ID e.g. 'CHEMBL1849'"
                        },
                        "activity_type": {
                            "type": "string",
                            "description": "Type of bioactivity measurement",
                            "enum": ["IC50", "Ki", "Kd", "EC50"],
                            "default": "IC50"
                        },
                        "max_compounds": {
                            "type": "integer",
                            "description": "Maximum compounds to retrieve (1-1000)",
                            "default": 300,
                            "minimum": 1,
                            "maximum": 1000
                        }
                    },
                    "required": ["target_chembl_id"]
                }
            },
            {
                "name": "get_chembl_compound",
                "description": (
                    "Get full molecule record from ChEMBL including canonical SMILES, "
                    "molecular formula, drug status, synonyms, and clinical phase."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "molecule_chembl_id": {
                            "type": "string",
                            "description": "ChEMBL molecule ID e.g. 'CHEMBL25'"
                        }
                    },
                    "required": ["molecule_chembl_id"]
                }
            }
        ]

    # ── Rate limiting ─────────────────────────────────────────────────────────

    def _wait(self) -> None:
        elapsed = time.time() - ChEMBLTools._last_call
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)
        ChEMBLTools._last_call = time.time()

    # ── Tool implementations ──────────────────────────────────────────────────

    def search_target(self, target_name: str, organism: str = "",
                      max_results: int = 5) -> ToolResult:
        try:
            from chembl_webresource_client.new_client import new_client
        except ImportError:
            return self._err("search_chembl_target", "chembl-webresource-client not installed.")

        try:
            self._wait()
            target_api = new_client.target
            filters: dict = {"target_synonym__icontains": target_name}
            if organism:
                filters["organism__icontains"] = organism

            raw = list(
                target_api.filter(**filters)
                .only([
                    "target_chembl_id", "pref_name", "target_type",
                    "organism", "target_components"
                ])[:max_results]
            )

            if not raw:
                return self._err(
                    "search_chembl_target",
                    f"No ChEMBL targets found for '{target_name}'"
                    + (f" in organism '{organism}'" if organism else "") + "."
                )

            targets = []
            for r in raw:
                targets.append({
                    "target_chembl_id": r.get("target_chembl_id", ""),
                    "pref_name": r.get("pref_name", ""),
                    "target_type": r.get("target_type", ""),
                    "organism": r.get("organism", ""),
                    "component_count": len(r.get("target_components") or []),
                    "chembl_url": f"{self._CHEMBL_URL}/target_report_card/{r.get('target_chembl_id', '')}/",
                })

            confidence = min(0.95, 0.6 + len(targets) * 0.07)
            return self._ok(
                tool_name="search_chembl_target",
                data={"targets": targets, "count": len(targets), "query": target_name},
                citations=[self._cite(
                    source="ChEMBL database (EMBL-EBI)",
                    url=f"{self._CHEMBL_URL}/target_report_card/{targets[0]['target_chembl_id']}/",
                    doi="10.1093/nar/gkac1049"
                )],
                confidence=confidence,
                rationale=f"ChEMBL is a manually curated database. {len(targets)} target(s) returned from live API query."
            )

        except Exception as exc:
            return self._err("search_chembl_target", str(exc))

    def fetch_activities(self, target_chembl_id: str,
                         activity_type: str = "IC50",
                         max_compounds: int = 300) -> ToolResult:
        try:
            from chembl_webresource_client.new_client import new_client
        except ImportError:
            return self._err("fetch_chembl_activities", "chembl-webresource-client not installed.")

        try:
            self._wait()
            activity_api = new_client.activity
            raw = list(
                activity_api.filter(
                    target_chembl_id=target_chembl_id,
                    standard_type=activity_type,
                    standard_relation="=",
                    standard_units="nM",
                ).only([
                    "molecule_chembl_id", "canonical_smiles", "standard_value",
                    "standard_type", "assay_chembl_id", "pchembl_value",
                    "data_validity_comment",
                ])[:max_compounds]
            )

            records = []
            skipped = 0
            for r in raw:
                if not r.get("standard_value") or not r.get("canonical_smiles"):
                    skipped += 1
                    continue
                try:
                    ic50 = float(r["standard_value"])
                    if ic50 <= 0:
                        skipped += 1
                        continue
                except (ValueError, TypeError):
                    skipped += 1
                    continue

                records.append({
                    "chembl_id": r.get("molecule_chembl_id", ""),
                    "smiles": r.get("canonical_smiles", ""),
                    "ic50_nm": ic50,
                    "activity_type": activity_type,
                    "assay_id": r.get("assay_chembl_id", ""),
                    "pchembl": float(r["pchembl_value"]) if r.get("pchembl_value") else None,
                    "data_validity": r.get("data_validity_comment") or "valid",
                })

            if not records:
                return self._err(
                    "fetch_chembl_activities",
                    f"No valid {activity_type} records with SMILES found for {target_chembl_id}."
                )

            n = len(records)
            active = sum(1 for r in records if r["ic50_nm"] < 1000)
            confidence = min(0.95, 0.50 + min(n, 100) / 200)

            warnings = []
            if skipped > 0:
                warnings.append(f"{skipped} records skipped (missing SMILES or value).")
            if n < 30:
                warnings.append(f"Only {n} compounds retrieved — QSAR model will have limited power.")

            return self._ok(
                tool_name="fetch_chembl_activities",
                data={
                    "target_chembl_id": target_chembl_id,
                    "activity_type": activity_type,
                    "total_records": n,
                    "active_count": active,
                    "inactive_count": n - active,
                    "skipped": skipped,
                    "compounds": records,
                },
                citations=[self._cite(
                    source=f"ChEMBL bioactivity data for target {target_chembl_id}",
                    url=f"{self._CHEMBL_URL}/target_report_card/{target_chembl_id}/",
                    doi="10.1093/nar/gkac1049"
                )],
                confidence=confidence,
                rationale=f"{n} real measured {activity_type} values from ChEMBL live API. Confidence scales with dataset size.",
                warnings=warnings,
            )

        except Exception as exc:
            return self._err("fetch_chembl_activities", str(exc))

    def get_compound(self, molecule_chembl_id: str) -> ToolResult:
        try:
            from chembl_webresource_client.new_client import new_client
        except ImportError:
            return self._err("get_chembl_compound", "chembl-webresource-client not installed.")

        try:
            self._wait()
            mol_api = new_client.molecule
            r = mol_api.get(molecule_chembl_id)

            if not r:
                return self._err("get_chembl_compound", f"{molecule_chembl_id} not found in ChEMBL.")

            props = r.get("molecule_properties") or {}
            record = {
                "chembl_id": r.get("molecule_chembl_id", ""),
                "pref_name": r.get("pref_name") or "",
                "smiles": r.get("molecule_structures", {}).get("canonical_smiles", ""),
                "inchi_key": r.get("molecule_structures", {}).get("standard_inchi_key", ""),
                "molecular_formula": props.get("full_molformula", ""),
                "mw": float(props.get("full_mwt") or 0),
                "alogp": float(props.get("alogp") or 0),
                "hbd": int(props.get("hbd") or 0),
                "hba": int(props.get("hba") or 0),
                "psa": float(props.get("psa") or 0),
                "rotatable_bonds": int(props.get("rtb") or 0),
                "drug_status": r.get("max_phase") or 0,
                "synonyms": [s.get("molecule_synonym", "") for s in (r.get("molecule_synonyms") or [])[:5]],
                "chembl_url": f"{self._CHEMBL_URL}/compound_report_card/{molecule_chembl_id}/",
            }

            return self._ok(
                tool_name="get_chembl_compound",
                data=record,
                citations=[self._cite(
                    source=f"ChEMBL compound {molecule_chembl_id}",
                    url=record["chembl_url"],
                    doi="10.1093/nar/gkac1049"
                )],
                confidence=0.95,
                rationale="ChEMBL compound records are manually curated from primary literature.",
            )

        except Exception as exc:
            return self._err("get_chembl_compound", str(exc))
