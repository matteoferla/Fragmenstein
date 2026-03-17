"""Direct SmallWorld API client with TLS 1.2 enforced.

sw.docking.org runs OpenSSL 1.1.1k which breaks on TLS 1.3 handshakes.
This client uses a raw ssl.SSLContext pinned to TLS 1.2 for all connections,
avoiding the monkey-patching issues with requests/urllib3.

API docs: https://wiki.docking.org/index.php/How_to_use_SmallWorld_API
"""

import json
import logging
import operator
import re
import ssl
import time
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

log = logging.getLogger(__name__)

BASE_URL = "https://sw.docking.org"

# TLS 1.2 context — created once, reused for all connections
_ssl_ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
_ssl_ctx.maximum_version = ssl.TLSVersion.TLSv1_2
_ssl_ctx.load_default_certs()

_DEFAULT_SUBMIT = {
    "dist": 8, "tdn": 6, "rdn": 6, "rup": 2, "ldn": 2, "lup": 2,
    "maj": 6, "min": 6, "sub": 6, "sdist": 12, "tup": 6,
    "scores": "Atom Alignment,ECFP4,Daylight",
}

# Column spec for results endpoint
_RESULT_COLUMNS = [
    "alignment", "dist", "ecfp4", "daylight", "topodist", "mces",
    "tdn", "tup", "rdn", "rup", "ldn", "lup", "mut", "maj", "min", "hyb", "sub",
]

_SEARCH_RANGES = {
    "dist": "0-12", "topodist": "0-8",
    "tdn": "0-6", "tup": "0-6", "rdn": "0-6", "rup": "0-2",
    "ldn": "0-2", "lup": "0-2", "maj": "0-6", "min": "0-6",
    "hyb": "0-6", "sub": "0-6",
}


def _build_view_params(hlid: int, start: int, length: int, draw: int) -> dict:
    """Build the column-spec params for /search/view."""
    params: dict[str, Any] = {"hlid": hlid, "start": start, "length": length, "draw": draw}
    for i, name in enumerate(_RESULT_COLUMNS):
        prefix = f"columns[{i}]"
        params[f"{prefix}[data]"] = str(i)
        params[f"{prefix}[name]"] = name
        params[f"{prefix}[searchable]"] = "true"
        params[f"{prefix}[orderable]"] = "true" if i > 0 else "false"
        params[f"{prefix}[search][value]"] = _SEARCH_RANGES.get(name, "")
        params[f"{prefix}[search][regex]"] = "false"
    params["order[0][column]"] = "0"
    params["order[0][dir]"] = "asc"
    params["search[value]"] = ""
    params["search[regex]"] = "false"
    return params


def _get(url: str, timeout: int = 600) -> bytes:
    """HTTP GET with TLS 1.2."""
    req = Request(url)
    with urlopen(req, context=_ssl_ctx, timeout=timeout) as resp:
        return resp.read()


def _submit_query(smiles: str, db: str, dist: int) -> int:
    """Submit a SMILES query, return the hit list ID (hlid)."""
    params = {**_DEFAULT_SUBMIT, "smi": smiles, "db": db, "dist": dist}
    url = f"{BASE_URL}/search/submit?{urlencode(params)}"

    raw = _get(url).decode("utf-8")
    hlid = -1
    for line in raw.splitlines():
        line = line.strip()
        if not line.startswith("data:"):
            continue
        try:
            datum = json.loads(line.removeprefix("data:").strip())
        except json.JSONDecodeError:
            continue
        if "hlid" in datum:
            hlid = datum["hlid"]
    if hlid == -1:
        raise ValueError(f"No hlid returned for SMILES: {smiles[:60]}")
    return hlid


def _get_results(hlid: int, length: int) -> pd.DataFrame:
    """Fetch results for a given hit list ID."""
    params = _build_view_params(hlid, start=0, length=length, draw=10)
    url = f"{BASE_URL}/search/view?{urlencode(params)}"

    data = json.loads(_get(url))
    if not data.get("recordsTotal"):
        return pd.DataFrame()

    rows = data.get("data", [])
    if not rows:
        return pd.DataFrame()

    # Column 0 is a nested dict (alignment info), rest are score columns
    df1 = pd.DataFrame(map(operator.itemgetter(0), rows))
    df2 = pd.DataFrame(rows).drop(columns=[0])
    df2.columns = _RESULT_COLUMNS[1:]
    df = pd.concat([df1, df2], axis=1)
    df["name"] = df.hitSmiles.str.split(expand=True)[1]
    df["smiles"] = df.hitSmiles.str.split(expand=True)[0]
    return df


def search_smiles(smiles: str, db: str, dist: int = 25, length: int = 200) -> pd.DataFrame:
    """Search SmallWorld for a single SMILES. Returns a DataFrame of hits."""
    hlid = _submit_query(smiles, db, dist)
    return _get_results(hlid, length)


def search_many(
    smiles_list: list[str],
    db: str,
    dist: int = 25,
    length: int = 200,
) -> pd.DataFrame:
    """Search SmallWorld for multiple SMILES with rate limiting."""
    results: list[pd.DataFrame] = []

    for i, smi in enumerate(smiles_list):
        log.info(f"SmallWorld: querying {i + 1}/{len(smiles_list)}: {smi[:60]}")
        try:
            df = search_smiles(smi, db=db, dist=dist, length=length)
            if not df.empty:
                df["query_index"] = i
                df["query_smiles"] = smi
                results.append(df)
                log.info(f"  → {len(df)} hits")
            else:
                log.info("  → 0 hits")
        except Exception as e:
            log.warning(f"  → failed: {e}")

        # Rate limit: sw.docking.org requires >= 5s between queries
        if i < len(smiles_list) - 1:
            time.sleep(5)

    if not results:
        return pd.DataFrame()
    return pd.concat(results, ignore_index=True)
