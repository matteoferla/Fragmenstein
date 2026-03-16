"""Wraps Laboratory.place() for the web API."""

import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem

from fragmenstein import Laboratory

from . import file_manager
from .combine_service import VICTOR_TYPES, load_session_hits, load_session_template
from .result_serializer import load_dataframe, save_dataframe

log = logging.getLogger(__name__)


def run_place(
    session_id: str,
    queries: pd.DataFrame | list[dict],
    victor_type: str = "Wictor",
    n_cores: int = -1,
    timeout: int = 240,
    covalent_resi: str | None = None,
) -> Path:
    """Run Laboratory.place() and save results. Returns path to results file."""
    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)

    # Resolve Victor type
    VictorClass = VICTOR_TYPES.get(victor_type)
    if VictorClass is None:
        from fragmenstein import Wictor
        VictorClass = Wictor

    # Configure
    Laboratory.Victor = VictorClass
    Laboratory.Victor.work_path = str(work_dir)
    Laboratory.Victor.monster_throw_on_discard = True
    Laboratory.Victor.error_to_catch = Exception

    lab = Laboratory(pdbblock=pdbblock, covalent_resi=covalent_resi)

    if isinstance(queries, list):
        queries = pd.DataFrame(queries)

    df = lab.place(queries, n_cores=n_cores, timeout=timeout)

    # Save results
    result_path = file_manager.results_dir(session_id, "place") / "results.pkl"
    save_dataframe(df, result_path)

    log.info(f"Place completed for session {session_id}: {len(df)} results")
    return result_path


def build_place_queries_from_similars(
    session_id: str,
    similars_df: pd.DataFrame,
    combine_result_path: str | None = None,
) -> pd.DataFrame:
    """Build placement queries from SmallWorld results + combine hits."""
    hits = load_session_hits(session_id)
    hit_map = {m.GetProp("_Name"): m for m in hits if m.HasProp("_Name")}

    queries = []
    for _, row in similars_df.iterrows():
        smiles = row.get("smiles", "")
        name = row.get("name", f"analog_{_}")

        # Map to original hits
        query_hits = hits  # Default to all hits
        if "hit_names" in row and isinstance(row["hit_names"], (list, tuple)):
            query_hits = [hit_map[n] for n in row["hit_names"] if n in hit_map]
            if not query_hits:
                query_hits = hits

        queries.append({"smiles": smiles, "name": name, "hits": query_hits})

    return pd.DataFrame(queries)
