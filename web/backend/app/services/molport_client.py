"""MolPort API client for similarity search."""

import logging
import threading

import requests

log = logging.getLogger(__name__)

# Max 3 concurrent chemical search queries
_CONCURRENCY_SEMAPHORE = threading.Semaphore(3)


class MolPortClient:
    """MolPort API v3.0 client for similarity search."""

    SEARCH_URL = "https://api.molport.com/api/chemical-search/search"

    def __init__(self, api_key: str):
        self._api_key = api_key

    def similarity_search(
        self,
        smiles: str,
        threshold: float = 0.8,
        max_results: int = 100,
        max_search_time: int = 60000,
    ) -> list[dict]:
        """Search MolPort for similar compounds.

        Returns list of dicts with keys: molport_id, smiles, similarity_index.
        """
        with _CONCURRENCY_SEMAPHORE:
            try:
                resp = requests.post(
                    self.SEARCH_URL,
                    json={
                        "API Key": self._api_key,
                        "Structure": smiles,
                        "Search Type": 4,  # similarity search
                        "Chemical Similarity Index": threshold,
                        "Maximum Search Time": max_search_time,
                        "Maximum Result Count": max_results,
                    },
                    timeout=max_search_time / 1000 + 10,
                )

                if resp.status_code == 429:
                    log.warning("MolPort: rate limited")
                    return []

                resp.raise_for_status()
                data = resp.json()

                status = data.get("Result", {}).get("Status", 0)
                if status != 1:
                    msg = data.get("Result", {}).get("Message", "Unknown error")
                    log.warning(f"MolPort: search failed for '{smiles[:40]}': {msg}")
                    return []

                molecules = data.get("Data", {}).get("Molecules", [])
                results = []
                for mol in molecules:
                    results.append({
                        "molport_id": mol.get("Molport Id", ""),
                        "smiles": mol.get("SMILES") or mol.get("Canonical SMILES", ""),
                        "similarity_index": mol.get("Similarity Index"),
                    })
                return results

            except requests.exceptions.Timeout:
                log.warning(f"MolPort: timeout for '{smiles[:40]}'")
                return []
            except Exception as e:
                log.warning(f"MolPort: error for '{smiles[:40]}': {e}")
                return []
