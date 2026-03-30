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
    use_originals: bool = True,
    combine_result_path: str | None = None,
) -> pd.DataFrame:
    """Build placement queries from SmallWorld/manual results + hits.

    use_originals=True: use original session hits for placement (default, recommended).
    use_originals=False: use the merger's unminimized mol as the hit template.
    """
    if similars_df.empty:
        return pd.DataFrame(columns=["smiles", "name", "hits"])

    hits = load_session_hits(session_id)
    if not hits and use_originals:
        raise ValueError("No hits loaded for session")

    # Load combine results for merger mols if use_originals=False
    combine_df = None
    if not use_originals and combine_result_path:
        combine_df = load_dataframe(Path(combine_result_path))

    # SmallWorld results may have hit_mols (set by similars_service) or not
    has_hit_mols = "hit_mols" in similars_df.columns

    queries = []
    seen = set()

    for _, row in similars_df.iterrows():
        smiles = row.get("smiles", "")
        if not smiles or not isinstance(smiles, str) or smiles in seen:
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        seen.add(smiles)
        name = row.get("name", f"analog_{len(queries)}")
        if not isinstance(name, str):
            name = str(name)

        if use_originals:
            # Use original hits
            query_hits = hits
            if has_hit_mols:
                row_hits = row.get("hit_mols")
                if isinstance(row_hits, (list, tuple)) and len(row_hits) > 0:
                    valid = [h for h in row_hits if h is not None and isinstance(h, Chem.Mol)]
                    if valid:
                        query_hits = valid
        else:
            # Use merger's unminimized mol as the hit template
            query_hits = hits  # fallback
            if combine_df is not None and "unminimized_mol" in combine_df.columns:
                # Try to find the merger that this analog came from
                query_smiles = row.get("query_smiles") or row.get("hitSmiles")
                if query_smiles and query_smiles in combine_df.get("simple_smiles", pd.Series()).values:
                    match = combine_df[combine_df["simple_smiles"] == query_smiles].iloc[0]
                    merger_mol = match.get("unminimized_mol")
                    if merger_mol is not None and isinstance(merger_mol, Chem.Mol):
                        query_hits = [merger_mol]

        queries.append({"smiles": smiles, "name": name, "hits": query_hits})

    if not queries:
        return pd.DataFrame(columns=["smiles", "name", "hits"])

    log.info(f"Built {len(queries)} placement queries (use_originals={use_originals})")
    return pd.DataFrame(queries)
