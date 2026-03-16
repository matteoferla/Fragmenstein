"""SmallWorld API wrapper for finding purchasable analogs."""

import logging
from pathlib import Path

import pandas as pd

from . import file_manager
from .result_serializer import load_dataframe, save_dataframe

log = logging.getLogger(__name__)


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
    from smallworld_api import SmallWorld

    # Load combine results
    combine_df = load_dataframe(Path(combine_result_path))

    # Filter to best results
    if outcome_filter:
        filtered = combine_df[combine_df["outcome"] == outcome_filter]
    else:
        filtered = combine_df[combine_df["outcome"].isin(["acceptable", "deviant"])]

    filtered = filtered.sort_values("∆∆G").head(top_n)

    if filtered.empty:
        log.warning("No results passed filter for SmallWorld search")
        result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
        save_dataframe(pd.DataFrame(), result_path)
        return result_path

    # Get SMILES to search
    smiles_col = "simple_smiles" if "simple_smiles" in filtered.columns else "smiles"
    smiles_list = filtered[smiles_col].dropna().tolist()

    # Run SmallWorld search
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
