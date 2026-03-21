"""PubMed/EuropePMC literature search — no API key required."""
from __future__ import annotations

import time
import httpx
from .base import BaseTool, ToolResult


class PubMedTools(BaseTool):

    _BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    _last_call: float = 0.0
    _delay: float = 0.4

    @classmethod
    def get_tool_definitions(cls) -> list[dict]:
        return [
            {
                "name": "search_literature",
                "description": (
                    "Search EuropePMC (PubMed-indexed) for peer-reviewed publications. "
                    "Returns titles, authors, journal, year, DOIs, and abstracts. "
                    "Use to find published evidence for targets, mechanisms, and compounds."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "Search terms e.g. 'InhA inhibitor tuberculosis' or 'EGFR kinase lung cancer'"
                        },
                        "max_results": {
                            "type": "integer",
                            "default": 10,
                            "minimum": 1,
                            "maximum": 50
                        },
                        "from_year": {
                            "type": "integer",
                            "default": 2010,
                            "description": "Filter to papers published from this year onward"
                        },
                        "open_access_only": {
                            "type": "boolean",
                            "default": False,
                            "description": "If true, return only open-access papers with full text"
                        }
                    },
                    "required": ["query"]
                }
            },
            {
                "name": "fetch_paper_details",
                "description": (
                    "Fetch full abstract and metadata for a specific paper by PubMed ID (PMID) or DOI. "
                    "Returns title, abstract, authors, journal, and citation metadata."
                ),
                "input_schema": {
                    "type": "object",
                    "properties": {
                        "pmid": {
                            "type": "string",
                            "description": "PubMed ID (numeric) or DOI string e.g. '10.1371/journal.pbio.3001721'"
                        }
                    },
                    "required": ["pmid"]
                }
            }
        ]

    def _wait(self) -> None:
        elapsed = time.time() - PubMedTools._last_call
        if elapsed < self._delay:
            time.sleep(self._delay - elapsed)
        PubMedTools._last_call = time.time()

    def _get(self, url: str, params: dict | None = None) -> dict | None:
        self._wait()
        try:
            r = httpx.get(url, params=params, timeout=15, follow_redirects=True)
            r.raise_for_status()
            return r.json()
        except Exception:
            raise

    def search_literature(self, query: str, max_results: int = 10,
                          from_year: int = 2010,
                          open_access_only: bool = False) -> ToolResult:
        try:
            q = f"({query}) AND FIRST_PDATE:[{from_year} TO 9999]"
            if open_access_only:
                q += " AND (OPEN_ACCESS:y)"

            params = {
                "query": q,
                "resulttype": "core",
                "format": "json",
                "pageSize": min(max_results, 50),
                "sort": "CITED desc",
            }
            data = self._get(f"{self._BASE}/search", params=params)

            if not data or not data.get("resultList", {}).get("result"):
                return self._err("search_literature",
                                 f"No EuropePMC results for query: '{query}'")

            papers = []
            for art in data["resultList"]["result"]:
                pmid = art.get("pmid") or art.get("id", "")
                doi = art.get("doi", "")
                papers.append({
                    "pmid": str(pmid),
                    "doi": doi,
                    "title": art.get("title", "").strip(),
                    "journal": art.get("journalTitle", ""),
                    "year": art.get("pubYear", ""),
                    "authors": art.get("authorString", ""),
                    "abstract": (art.get("abstractText") or "")[:500],
                    "citation_count": art.get("citedByCount", 0),
                    "is_open_access": art.get("isOpenAccess", "N") == "Y",
                    "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                    "europepmc_url": f"https://europepmc.org/article/MED/{pmid}" if pmid else "",
                })

            n = len(papers)
            confidence = 0.85 if n >= 5 else 0.60

            return self._ok(
                tool_name="search_literature",
                data={"query": query, "total_found": data.get("hitCount", n),
                      "returned": n, "papers": papers},
                citations=[self._cite(
                    source="EuropePMC / PubMed",
                    url=f"https://europepmc.org/search?query={query.replace(' ', '+')}",
                    doi=""
                )],
                confidence=confidence,
                rationale=f"{n} peer-reviewed papers found; EuropePMC indexes PubMed. Citation count used for relevance ranking.",
            )

        except Exception as exc:
            return self._err("search_literature", str(exc))

    def fetch_paper_details(self, pmid: str) -> ToolResult:
        try:
            # Strip DOI prefix if present
            id_val = pmid.replace("https://doi.org/", "").strip()
            if id_val.startswith("10."):
                source = "DOI"
                params = {"query": f"DOI:{id_val}", "resulttype": "core",
                          "format": "json", "pageSize": 1}
                data = self._get(f"{self._BASE}/search", params=params)
                if not data or not data.get("resultList", {}).get("result"):
                    return self._err("fetch_paper_details", f"No paper found for DOI: {id_val}")
                art = data["resultList"]["result"][0]
            else:
                params = {"query": f"EXT_ID:{id_val} AND SRC:MED", "resulttype": "core",
                          "format": "json", "pageSize": 1}
                data = self._get(f"{self._BASE}/search", params=params)
                if not data or not data.get("resultList", {}).get("result"):
                    return self._err("fetch_paper_details", f"No paper found for PMID: {id_val}")
                art = data["resultList"]["result"][0]

            paper_pmid = art.get("pmid") or art.get("id", "")
            doi = art.get("doi", "")
            record = {
                "pmid": str(paper_pmid),
                "doi": doi,
                "title": art.get("title", "").strip(),
                "journal": art.get("journalTitle", ""),
                "year": art.get("pubYear", ""),
                "authors": art.get("authorString", ""),
                "abstract": art.get("abstractText") or "Abstract not available.",
                "keywords": art.get("keywordList", {}).get("keyword", []),
                "citation_count": art.get("citedByCount", 0),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{paper_pmid}/",
            }

            return self._ok(
                tool_name="fetch_paper_details",
                data=record,
                citations=[self._cite(
                    source=f"{record['authors'].split(',')[0]} et al., {record['journal']} {record['year']}",
                    url=record["pubmed_url"],
                    doi=doi
                )],
                confidence=0.90,
                rationale="PubMed-indexed article; full metadata retrieved from EuropePMC.",
            )

        except Exception as exc:
            return self._err("fetch_paper_details", str(exc))
