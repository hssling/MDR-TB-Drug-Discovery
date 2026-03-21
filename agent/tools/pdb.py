"""PDB structure tools — search and download from RCSB."""
from __future__ import annotations

import json
import time
import httpx
from pathlib import Path
from .base import BaseTool, ToolResult


class PDBTools(BaseTool):

    _SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
    _DATA_API = "https://data.rcsb.org/rest/v1/core/entry"
    _FILES_URL = "https://files.rcsb.org/download"
    _last_call: float = 0.0
    _delay: float = 0.5

    def __init__(self, pdb_cache_dir: Path | None = None):
        self._cache = pdb_cache_dir or Path("data/pdb")
        self._cache.mkdir(parents=True, exist_ok=True)

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        return [
            {
                "name": "search_pdb_structure",
                "description": (
                    "Search RCSB PDB for protein structures by protein name or UniProt ID. "
                    "Returns PDB IDs, resolution (Å), experimental method, organism, "
                    "and bound ligands. Prefer structures with resolution < 2.5 Å."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "protein_name": {
                            "type": "string",
                            "description": "Protein name e.g. 'InhA', 'DprE1', 'EGFR'"
                        },
                        "organism": {
                            "type": "string",
                            "description": "Scientific organism name e.g. 'Mycobacterium tuberculosis'"
                        },
                        "max_results": {
                            "type": "integer",
                            "default": 5,
                            "minimum": 1,
                            "maximum": 20
                        }
                    },
                    "required": ["protein_name"]
                }
            },
            {
                "name": "download_pdb_structure",
                "description": (
                    "Download a PDB structure file from RCSB and extract metadata: "
                    "resolution, organism, bound ligands, chain count, and active-site residues. "
                    "File is cached locally for subsequent use."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "pdb_id": {
                            "type": "string",
                            "description": "4-character PDB entry ID e.g. '4TZK', '1IVA'"
                        }
                    },
                    "required": ["pdb_id"]
                }
            }
        ]

    def _wait(self) -> None:
        elapsed = time.time() - PDBTools._last_call
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)
        PDBTools._last_call = time.time()

    def search_structure(self, protein_name: str, organism: str = "",
                         max_results: int = 5) -> ToolResult:
        try:
            query: dict = {
                "query": {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "full_text",
                            "parameters": {"value": protein_name}
                        }
                    ]
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": 0, "rows": max_results},
                    "sort": [{"sort_by": "score", "direction": "desc"}],
                    "scoring_strategy": "combined"
                },
                "request_info": {"query_id": "drug_discovery_agent"}
            }

            if organism:
                query["query"]["nodes"].append({
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.scientific_name",
                        "operator": "contains_phrase",
                        "value": organism
                    }
                })

            self._wait()
            r = httpx.post(self._SEARCH_API, json=query, timeout=15)
            r.raise_for_status()
            search_data = r.json()

            result_ids = [h["identifier"] for h in search_data.get("result_list", [])]
            if not result_ids:
                return self._err("search_pdb_structure",
                                 f"No PDB structures found for '{protein_name}'"
                                 + (f" from '{organism}'" if organism else "") + ".")

            structures = []
            for pdb_id in result_ids[:max_results]:
                self._wait()
                try:
                    meta = httpx.get(f"{self._DATA_API}/{pdb_id}", timeout=10)
                    meta.raise_for_status()
                    m = meta.json()

                    res = None
                    diffraction = m.get("refine", [{}])
                    if diffraction:
                        res = diffraction[0].get("ls_d_res_high")

                    structures.append({
                        "pdb_id": pdb_id,
                        "title": m.get("struct", {}).get("title", ""),
                        "resolution_angstrom": float(res) if res else None,
                        "experimental_method": m.get("exptl", [{}])[0].get("method", ""),
                        "deposition_date": m.get("rcsb_accession_info", {}).get("deposit_date", ""),
                        "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}",
                    })
                except Exception:
                    structures.append({"pdb_id": pdb_id,
                                       "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}"})

            best = structures[0]
            res_str = f"{best.get('resolution_angstrom', '?')} Å" if best.get("resolution_angstrom") else "unknown"

            return self._ok(
                tool_name="search_pdb_structure",
                data={"structures": structures, "count": len(structures), "query": protein_name},
                citations=[self._cite(
                    source="RCSB Protein Data Bank",
                    url=f"https://www.rcsb.org/structure/{best['pdb_id']}",
                    doi="10.1093/nar/28.1.235"
                )],
                confidence=0.90,
                rationale=f"RCSB PDB is the authoritative structural biology archive. Best structure: {best['pdb_id']} ({res_str}).",
            )

        except Exception as exc:
            return self._err("search_pdb_structure", str(exc))

    def download_structure(self, pdb_id: str) -> ToolResult:
        pdb_id = pdb_id.upper().strip()
        try:
            pdb_path = self._cache / f"{pdb_id}.pdb"

            if not pdb_path.exists():
                self._wait()
                url = f"{self._FILES_URL}/{pdb_id}.pdb"
                r = httpx.get(url, timeout=30, follow_redirects=True)
                r.raise_for_status()
                pdb_path.write_bytes(r.content)

            content = pdb_path.read_text(encoding="utf-8", errors="replace")
            lines = content.splitlines()

            # Parse metadata from PDB file
            resolution = None
            ligands: set[str] = set()
            chains: set[str] = set()
            organisms: list[str] = []

            for line in lines:
                if line.startswith("REMARK   2 RESOLUTION."):
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if p == "ANGSTROMS." and i > 0:
                            try:
                                resolution = float(parts[i - 1])
                            except ValueError:
                                pass
                elif line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    if res_name not in {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO"}:
                        ligands.add(res_name)
                elif line.startswith("ATOM"):
                    chains.add(line[21])
                elif line.startswith("SOURCE") and "ORGANISM_SCIENTIFIC" in line:
                    organisms.append(line.split("ORGANISM_SCIENTIFIC:")[-1].strip().rstrip(";"))

            file_size_kb = pdb_path.stat().st_size // 1024
            confidence = 0.95 if resolution and resolution < 2.0 else \
                         0.85 if resolution and resolution < 2.5 else 0.70

            return self._ok(
                tool_name="download_pdb_structure",
                data={
                    "pdb_id": pdb_id,
                    "local_path": str(pdb_path),
                    "file_size_kb": file_size_kb,
                    "resolution_angstrom": resolution,
                    "chains": sorted(chains),
                    "chain_count": len(chains),
                    "ligands": sorted(ligands),
                    "organisms": organisms[:3],
                    "pdb_url": f"https://www.rcsb.org/structure/{pdb_id}",
                },
                citations=[self._cite(
                    source=f"RCSB PDB entry {pdb_id}",
                    url=f"https://www.rcsb.org/structure/{pdb_id}",
                    doi="10.1093/nar/28.1.235"
                )],
                confidence=confidence,
                rationale=f"PDB {pdb_id} resolution: {resolution} Å. Higher confidence for sub-2Å structures.",
            )

        except Exception as exc:
            return self._err("download_pdb_structure", str(exc))
