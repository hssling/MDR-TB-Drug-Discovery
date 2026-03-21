
"""PubChem tools — compound search and bioassay retrieval via PUG REST API."""
from __future__ import annotations

import time
import json
import httpx
from .base import BaseTool, ToolResult


class PubChemTools(BaseTool):

    _BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    _last_call: float = 0.0
    _delay: float = 0.3

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        return [
            {
                "name": "search_pubchem_compound",
                "description": (
                    "Search PubChem for compounds by name or CID. Returns canonical SMILES, "
                    "molecular formula, MW, CID, and IUPAC name. Useful for looking up "
                    "reference compounds, approved drugs, or known bioactive molecules."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "Compound name (e.g. 'aspirin', 'triclosan', 'bedaquiline') or PubChem CID"
                        },
                        "query_type": {
                            "type": "string",
                            "enum": ["name", "cid"],
                            "default": "name",
                            "description": "Whether query is a compound name or numeric CID"
                        },
                        "max_results": {
                            "type": "integer",
                            "default": 5,
                            "minimum": 1,
                            "maximum": 20
                        }
                    },
                    "required": ["query"]
                }
            },
            {
                "name": "get_pubchem_bioassays",
                "description": (
                    "Retrieve PubChem BioAssay activity summary for a compound CID. "
                    "Returns active/inactive assay counts, tested concentration ranges, "
                    "and activity outcome summaries."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "cid": {
                            "type": "integer",
                            "description": "PubChem compound CID (numeric identifier)"
                        },
                        "max_assays": {
                            "type": "integer",
                            "default": 20,
                            "minimum": 1,
                            "maximum": 100
                        }
                    },
                    "required": ["cid"]
                }
            }
        ]

    def _wait(self) -> None:
        elapsed = time.time() - PubChemTools._last_call
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)
        PubChemTools._last_call = time.time()

    def _get(self, url: str, timeout: int = 15) -> dict | list | None:
        self._wait()
        try:
            r = httpx.get(url, timeout=timeout, follow_redirects=True)
            r.raise_for_status()
            return r.json()
        except httpx.HTTPStatusError as e:
            if e.response.status_code == 404:
                return None
            raise
        except Exception:
            raise

    def search_compound(self, query: str, query_type: str = "name",
                        max_results: int = 5) -> ToolResult:
        props = "CanonicalSMILES,MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES"
        try:
            if query_type == "cid":
                url = f"{self._BASE}/compound/cid/{query}/property/{props}/JSON"
            else:
                q = query.replace(" ", "%20")
                url = f"{self._BASE}/compound/name/{q}/property/{props}/JSON"

            data = self._get(url)
            if not data:
                return self._err("search_pubchem_compound", f"No PubChem compound found for '{query}'.")

            compounds = []
            table = data.get("PropertyTable", {}).get("Properties", [])
            for row in table[:max_results]:
                cid = row.get("CID")
                compounds.append({
                    "cid": cid,
                    "canonical_smiles": row.get("CanonicalSMILES", ""),
                    "isomeric_smiles": row.get("IsomericSMILES", ""),
                    "molecular_formula": row.get("MolecularFormula", ""),
                    "mw": float(row.get("MolecularWeight") or 0),
                    "iupac_name": row.get("IUPACName", ""),
                    "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                })

            if not compounds:
                return self._err("search_pubchem_compound", f"Empty result for '{query}'.")

            return self._ok(
                tool_name="search_pubchem_compound",
                data={"query": query, "compounds": compounds, "count": len(compounds)},
                citations=[self._cite(
                    source="PubChem (NCBI)",
                    url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{compounds[0]['cid']}",
                    doi="10.1093/nar/gkac956"
                )],
                confidence=0.90,
                rationale="PubChem data is curated by NCBI from depositor submissions and literature.",
            )

        except Exception as exc:
            return self._err("search_pubchem_compound", str(exc))

    def get_bioassays(self, cid: int, max_assays: int = 20) -> ToolResult:
        try:
            url = f"{self._BASE}/compound/cid/{cid}/assaysummary/JSON"
            data = self._get(url)
            if not data:
                return self._err("get_pubchem_bioassays", f"No bioassay data found for CID {cid}.")

            table = data.get("Table", {})
            cols = table.get("Columns", {}).get("Column", [])
            rows = table.get("Row", [])

            assays = []
            for row in rows[:max_assays]:
                cells = row.get("Cell", [])
                record = dict(zip(cols, cells))
                assays.append({
                    "assay_id": record.get("AID", ""),
                    "outcome": record.get("Activity Outcome", ""),
                    "target_name": record.get("Target Name", ""),
                    "assay_name": record.get("Assay Name", ""),
                })

            active = sum(1 for a in assays if str(a.get("outcome", "")).lower() == "active")

            return self._ok(
                tool_name="get_pubchem_bioassays",
                data={
                    "cid": cid,
                    "total_assays_shown": len(assays),
                    "active_count": active,
                    "assays": assays,
                    "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Biological-Test-Results",
                },
                citations=[self._cite(
                    source=f"PubChem BioAssay data for CID {cid}",
                    url=f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
                    doi="10.1093/nar/gkac956"
                )],
                confidence=0.80,
                rationale="PubChem BioAssay data is depositor-submitted; curated but not uniformly peer-reviewed.",
            )

        except Exception as exc:
            return self._err("get_pubchem_bioassays", str(exc))
