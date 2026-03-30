"""ChemSpace API client with token management."""

import logging
import threading
import time

import requests

log = logging.getLogger(__name__)

# Rate limit: 40 req/min → 1.5s between requests
_RATE_LIMIT_DELAY = 1.5
_TOKEN_REFRESH_MARGIN = 30  # seconds before expiry to refresh


class ChemSpaceClient:
    """Thread-safe ChemSpace API client with automatic token refresh."""

    BASE_URL = "https://api.chem-space.com"

    def __init__(self, api_key: str):
        self._api_key = api_key
        self._token: str | None = None
        self._token_expires_at: float = 0
        self._lock = threading.Lock()
        self._last_request_time: float = 0

    def _ensure_token(self):
        """Get or refresh the access token (thread-safe)."""
        with self._lock:
            if self._token and time.time() < self._token_expires_at - _TOKEN_REFRESH_MARGIN:
                return
            resp = requests.get(
                f"{self.BASE_URL}/auth/token",
                headers={
                    "Authorization": f"Bearer {self._api_key}",
                    "Accept": "application/json",
                },
                timeout=15,
            )
            resp.raise_for_status()
            data = resp.json()
            self._token = data["access_token"]
            self._token_expires_at = time.time() + data.get("expires_in", 3600)
            log.info("ChemSpace: access token refreshed")

    def _rate_limit(self):
        """Enforce rate limiting between requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < _RATE_LIMIT_DELAY:
            time.sleep(_RATE_LIMIT_DELAY - elapsed)
        self._last_request_time = time.time()

    def similarity_search(
        self,
        smiles: str,
        count: int = 50,
        categories: str = "CSSS,CSMS",
        ship_to_country: str = "US",
    ) -> list[dict]:
        """Search ChemSpace for similar compounds.

        Returns list of dicts with keys: csId, smiles, molecular_weight, logP, TPSA.
        """
        self._ensure_token()
        self._rate_limit()

        try:
            resp = requests.post(
                f"{self.BASE_URL}/v4/search/sim",
                headers={
                    "Authorization": f"Bearer {self._token}",
                    "Accept": "application/json; version=4.1",
                },
                params={
                    "shipToCountry": ship_to_country,
                    "count": count,
                    "page": 1,
                    "categories": categories,
                },
                files={"SMILES": (None, smiles)},
                timeout=30,
            )

            if resp.status_code == 429:
                log.warning("ChemSpace: rate limited, waiting 5s")
                time.sleep(5)
                return self.similarity_search(smiles, count, categories, ship_to_country)

            if resp.status_code == 401:
                # Token expired, refresh and retry once
                self._token = None
                self._ensure_token()
                return self.similarity_search(smiles, count, categories, ship_to_country)

            if resp.status_code == 400:
                log.warning(f"ChemSpace: invalid SMILES '{smiles[:40]}', skipping")
                return []

            resp.raise_for_status()
            data = resp.json()
            items = data.get("items", [])

            results = []
            for item in items:
                results.append({
                    "csId": item.get("csId", ""),
                    "smiles": item.get("SMILES", ""),
                    "molecular_weight": item.get("MW"),
                    "logP": item.get("logP"),
                    "TPSA": item.get("TPSA"),
                })
            return results

        except requests.exceptions.Timeout:
            log.warning(f"ChemSpace: timeout for '{smiles[:40]}'")
            return []
        except Exception as e:
            log.warning(f"ChemSpace: error for '{smiles[:40]}': {e}")
            return []
