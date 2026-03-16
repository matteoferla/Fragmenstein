"""Wraps Laboratory.place() for the web API."""

import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem

from . import file_manager
from .combine_service import VICTOR_TYPES, WebLaboratory, load_session_hits, load_session_template
from .result_serializer import load_dataframe, save_dataframe

log = logging.getLogger(__name__)


def run_place(
    session_id: str,
    queries: pd.DataFrame | list[dict],
    victor_type: str = "Wictor",
    n_cores: int = -1,
    timeout: int = 240,
    merging_mode: str = "expansion",
    run_plip: bool = False,
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
    WebLaboratory.Victor = VictorClass
    WebLaboratory.Victor.work_path = str(work_dir)
    WebLaboratory.Victor.monster_throw_on_discard = True
    WebLaboratory.Victor.error_to_catch = Exception

    lab = WebLaboratory(pdbblock=pdbblock, covalent_resi=covalent_resi, run_plip=run_plip)
    lab._merging_mode = merging_mode

    if isinstance(queries, list):
        queries = pd.DataFrame(queries)

    if queries.empty:
        raise ValueError("No placement queries — the SmallWorld search returned no analogs. "
                         "Try relaxing the outcome filter or increasing the distance parameter.")

    df = lab.place(queries, n_cores=n_cores, timeout=timeout)

    # Save results
    result_path = file_manager.results_dir(session_id, "place") / "results.pkl"
    save_dataframe(df, result_path)

    log.info(f"Place completed for session {session_id}: {len(df)} results")
    return result_path


def build_place_queries_from_similars(
    session_id: str,
    similars_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build placement queries from SmallWorld results + session hits.

    SmallWorld results have columns: smiles, name, hitSmiles, topodist, etc.
    We need to build: smiles, name, hits (List[Mol]) for Laboratory.place().
    """
    if similars_df.empty:
        return pd.DataFrame(columns=["smiles", "name", "hits"])

    hits = load_session_hits(session_id)
    if not hits:
        raise ValueError("No hits loaded for session")

    # SmallWorld results may have hit_mols (set by similars_service) or not
    has_hit_mols = "hit_mols" in similars_df.columns

    queries = []
    seen = set()  # Deduplicate by SMILES

    for _, row in similars_df.iterrows():
        smiles = row.get("smiles", "")
        if not smiles or not isinstance(smiles, str) or smiles in seen:
            continue

        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        seen.add(smiles)
        name = row.get("name", f"analog_{len(queries)}")
        if not isinstance(name, str):
            name = str(name)

        # Determine which hits to use for this analog
        query_hits = hits  # default: all session hits
        if has_hit_mols:
            row_hits = row.get("hit_mols")
            if isinstance(row_hits, (list, tuple)) and len(row_hits) > 0:
                valid = [h for h in row_hits if h is not None and isinstance(h, Chem.Mol)]
                if valid:
                    query_hits = valid

        queries.append({"smiles": smiles, "name": name, "hits": query_hits})

    if not queries:
        return pd.DataFrame(columns=["smiles", "name", "hits"])

    log.info(f"Built {len(queries)} placement queries from {len(similars_df)} SmallWorld results")
    return pd.DataFrame(queries)
