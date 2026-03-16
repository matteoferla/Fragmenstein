"""SmallWorld API wrapper for finding purchasable analogs.

Root cause of SSL errors: sw.docking.org runs Apache/OpenSSL 1.1.1k which has a
broken TLS 1.3 handshake. Anaconda's OpenSSL 3.5.0 negotiates TLS 1.3 by default,
causing UNEXPECTED_EOF_WHILE_READING. Fix: cap ssl.SSLContext to TLS 1.2 for the
duration of the SmallWorld call.
"""

import logging
import ssl
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

from . import file_manager
from .result_serializer import load_dataframe, save_dataframe

log = logging.getLogger(__name__)


@contextmanager
def _tls12_only():
    """Temporarily cap all new SSLContext instances to TLS 1.2.

    sw.docking.org's OpenSSL 1.1.1k breaks on TLS 1.3 handshakes with newer
    clients. This context manager patches ssl.create_default_context so any
    SSL context created (by requests/urllib3/smallworld_api) is capped at
    TLS 1.2, then restores the original after the block completes.
    """
    # Patch create_default_context (used by urllib3 internally)
    original_create = ssl.create_default_context

    def patched_create(*args, **kwargs):
        ctx = original_create(*args, **kwargs)
        ctx.maximum_version = ssl.TLSVersion.TLSv1_2
        return ctx

    ssl.create_default_context = patched_create

    # Also patch urllib3's context creator directly
    try:
        from urllib3.util.ssl_ import create_urllib3_context as _orig_u3
        import urllib3.util.ssl_

        def patched_u3(*args, **kwargs):
            ctx = _orig_u3(*args, **kwargs)
            ctx.maximum_version = ssl.TLSVersion.TLSv1_2
            return ctx

        urllib3.util.ssl_.create_urllib3_context = patched_u3
        _has_u3 = True
    except ImportError:
        _has_u3 = False

    try:
        yield
    finally:
        ssl.create_default_context = original_create
        if _has_u3:
            urllib3.util.ssl_.create_urllib3_context = _orig_u3


def run_similars_search(
    session_id: str,
    combine_result_path: str,
    top_n: int = 100,
    dist: int = 25,
    length: int = 200,
    db: str = "REAL_dataset",
    outcome_filter: str = "acceptable",
) -> Path:
    """Search SmallWorld for purchasable analogs of combine results."""
    # Load combine results
    combine_df = load_dataframe(Path(combine_result_path))

    # Filter to best results
    if outcome_filter:
        filtered = combine_df[combine_df["outcome"] == outcome_filter]
    else:
        # Empty filter = all outcomes
        filtered = combine_df.copy()

    # Sort by energy (need valid ddG values on top)
    if "∆∆G" in filtered.columns:
        filtered = filtered.sort_values("∆∆G", na_position="last").head(top_n)
    else:
        filtered = filtered.head(top_n)

    if filtered.empty:
        log.warning("No results passed filter for SmallWorld search")
        result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
        save_dataframe(pd.DataFrame(), result_path)
        return result_path

    # Get SMILES to search
    smiles_col = "simple_smiles" if "simple_smiles" in filtered.columns else "smiles"
    smiles_list = filtered[smiles_col].dropna().tolist()

    # Run SmallWorld search with TLS 1.2 (sw.docking.org breaks on TLS 1.3)
    with _tls12_only():
        from smallworld_api import SmallWorld

        sws = SmallWorld()
        sw_db = getattr(sws, db, sws.REAL_dataset)

        results = sws.search_many(
            smiles_list,
            dist=dist,
            length=length,
            db=sw_db,
            tolerated_exceptions=Exception,
        )

    # Map back to original hit molecules
    if "hit_mols" in filtered.columns:
        hit_mols_map = dict(zip(filtered[smiles_col], filtered["hit_mols"]))
        if not results.empty and "hitSmiles" in results.columns:
            results["hits"] = results["hitSmiles"].map(
                lambda s: hit_mols_map.get(s, [])
            )

    # Save results
    result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
    save_dataframe(results, result_path)

    log.info(f"SmallWorld search completed: {len(results)} analogs found")
    return result_path
