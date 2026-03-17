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


def parse_manual_smiles(session_id: str, smiles_text: str) -> tuple[Path, int, int]:
    """Parse pasted SMILES text into a similars DataFrame.

    Format: one compound per line, SMILES<whitespace>NAME.
    Returns (result_path, valid_count, invalid_count).
    """
    from rdkit import Chem

    rows = []
    invalid = 0

    for i, line in enumerate(smiles_text.strip().split("\n")):
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)  # split on any whitespace, max 2 parts
        smiles = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"compound_{i + 1}"

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid += 1
            continue

        rows.append({"smiles": Chem.MolToSmiles(mol), "name": name})

    df = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["smiles", "name"])
    result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
    save_dataframe(df, result_path)
    return result_path, len(rows), invalid


def parse_uploaded_file(session_id: str, filepath: Path) -> tuple[Path, int, int]:
    """Parse a CSV/TSV/Excel file into a similars DataFrame.

    Auto-detects SMILES and name columns. Preserves extra columns.
    Returns (result_path, valid_count, invalid_count).
    """
    from rdkit import Chem

    suffix = filepath.suffix.lower()
    if suffix in (".csv", ".tsv"):
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(filepath, sep=sep)
    elif suffix == ".xlsx":
        df = pd.read_excel(filepath)
    else:
        raise ValueError(f"Unsupported file format: {suffix}. Use .csv, .tsv, or .xlsx")

    # Auto-detect SMILES column
    smiles_col = None
    for candidate in ["smiles", "SMILES", "Smiles", "canonical_smiles", "CanonicalSMILES", "smi"]:
        if candidate in df.columns:
            smiles_col = candidate
            break
    if smiles_col is None:
        raise ValueError(f"No SMILES column found. Expected one of: smiles, SMILES, canonical_smiles. Got: {list(df.columns)}")

    # Auto-detect name column
    name_col = None
    for candidate in ["name", "Name", "NAME", "id", "ID", "compound_id", "vendor_id", "Vendor_ID", "title"]:
        if candidate in df.columns:
            name_col = candidate
            break

    # Validate SMILES and build clean DataFrame
    valid_rows = []
    invalid = 0

    for i, row in df.iterrows():
        raw_smiles = str(row[smiles_col]).strip()
        mol = Chem.MolFromSmiles(raw_smiles)
        if mol is None:
            invalid += 1
            continue

        new_row = dict(row)
        new_row["smiles"] = Chem.MolToSmiles(mol)
        if name_col:
            new_row["name"] = str(row[name_col])
        elif "name" not in new_row:
            new_row["name"] = f"compound_{i + 1}"
        valid_rows.append(new_row)

    result_df = pd.DataFrame(valid_rows) if valid_rows else pd.DataFrame(columns=["smiles", "name"])
    result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
    save_dataframe(result_df, result_path)
    return result_path, len(valid_rows), invalid


def filter_library_against_mergers(
    session_id: str,
    library_path: Path,
    combine_result_path: str,
    top_n: int = 200,
    outcome_filter: str = "acceptable",
) -> tuple[Path, int, int, int]:
    """Filter a large compound library by Tanimoto similarity to combine mergers.

    Uses Morgan fingerprints (radius=2) for fast similarity calculation.
    Returns (result_path, library_size, filtered_count, merger_count).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs

    # Load library
    library_df = load_dataframe(library_path)
    if library_df.empty:
        raise ValueError("Library is empty")

    library_size = len(library_df)

    # Load combine results for merger SMILES
    combine_df = load_dataframe(Path(combine_result_path))
    if outcome_filter:
        combine_df = combine_df[combine_df["outcome"] == outcome_filter]
    else:
        combine_df = combine_df.copy()

    if combine_df.empty:
        raise ValueError("No mergers found matching the outcome filter")

    smiles_col = "simple_smiles" if "simple_smiles" in combine_df.columns else "smiles"
    merger_smiles = combine_df[smiles_col].dropna().tolist()

    # Build Morgan fingerprints for mergers
    merger_fps = []
    for smi in merger_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            merger_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
    if not merger_fps:
        raise ValueError("No valid merger fingerprints could be generated")

    log.info(f"Filtering {library_size} library compounds against {len(merger_fps)} mergers")

    # Compute max Tanimoto similarity for each library compound against all mergers
    scores = []
    smiles_key = None
    for candidate in ["smiles", "SMILES", "Smiles", "canonical_smiles"]:
        if candidate in library_df.columns:
            smiles_key = candidate
            break
    if smiles_key is None:
        raise ValueError("No SMILES column in library")

    for idx, row in library_df.iterrows():
        smi = str(row[smiles_key])
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            scores.append(0.0)
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        max_sim = max(DataStructs.TanimotoSimilarity(fp, mfp) for mfp in merger_fps)
        scores.append(max_sim)

    library_df["tanimoto_to_merger"] = scores

    # Sort by similarity, take top N
    filtered = library_df.sort_values("tanimoto_to_merger", ascending=False).head(top_n)

    # Ensure standard column names
    if smiles_key != "smiles":
        filtered = filtered.rename(columns={smiles_key: "smiles"})

    result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
    save_dataframe(filtered, result_path)

    log.info(f"Filtered {library_size} → {len(filtered)} compounds (top {top_n} by Tanimoto)")
    return result_path, library_size, len(filtered), len(merger_fps)


def run_pubchem_search(
    session_id: str,
    combine_result_path: str,
    top_n: int = 50,
    threshold: int = 80,
    max_per_query: int = 20,
    outcome_filter: str = "acceptable",
) -> Path:
    """Search PubChem for similar compounds using fingerprint similarity."""
    import time
    import requests
    from rdkit import Chem

    combine_df = load_dataframe(Path(combine_result_path))

    if outcome_filter:
        filtered = combine_df[combine_df["outcome"] == outcome_filter]
    else:
        filtered = combine_df.copy()

    if "∆∆G" in filtered.columns:
        filtered = filtered.sort_values("∆∆G", na_position="last").head(top_n)
    else:
        filtered = filtered.head(top_n)

    if filtered.empty:
        result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
        save_dataframe(pd.DataFrame(), result_path)
        return result_path

    smiles_col = "simple_smiles" if "simple_smiles" in filtered.columns else "smiles"
    smiles_list = filtered[smiles_col].dropna().tolist()

    all_results = []
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    log.info(f"PubChem: searching {len(smiles_list)} SMILES, threshold={threshold}")

    for i, smi in enumerate(smiles_list):
        if not smi or not smi.strip():
            continue
        # Clean SMILES: remove explicit H for PubChem compatibility
        clean_mol = Chem.MolFromSmiles(smi)
        if clean_mol is None:
            log.warning(f"PubChem: invalid query SMILES '{smi}', skipping")
            continue
        clean_smi = Chem.MolToSmiles(clean_mol)

        log.info(f"PubChem: querying {i + 1}/{len(smiles_list)}: {clean_smi[:50]}")
        try:
            url = (
                f"{base_url}/fastsimilarity_2d/smiles/{requests.utils.quote(clean_smi, safe='')}"
                f"/property/CanonicalSMILES,IUPACName,MolecularWeight/JSON"
                f"?Threshold={threshold}&MaxRecords={max_per_query}"
            )
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                props = data.get("PropertyTable", {}).get("Properties", [])
                for p in props:
                    canon = p.get("CanonicalSMILES") or p.get("ConnectivitySMILES") or p.get("IsomericSMILES", "")
                    if not canon:
                        continue
                    mol = Chem.MolFromSmiles(canon)
                    if mol is None:
                        continue
                    all_results.append({
                        "smiles": Chem.MolToSmiles(mol),
                        "name": f"CID-{p.get('CID', '?')}",
                        "molecular_weight": p.get("MolecularWeight"),
                        "iupac_name": p.get("IUPACName", ""),
                        "query_smiles": clean_smi,
                    })
                log.info(f"  → {len(props)} hits")
            elif resp.status_code == 404:
                log.info(f"  → 0 hits")
            else:
                log.warning(f"PubChem returned {resp.status_code} for '{clean_smi[:40]}'")
        except Exception as e:
            log.warning(f"PubChem search failed for '{clean_smi[:40]}': {e}")

        # Rate limit: PubChem allows 5 requests/second
        time.sleep(0.25)

    results = pd.DataFrame(all_results) if all_results else pd.DataFrame()
    result_path = file_manager.results_dir(session_id, "similars") / "results.pkl"
    save_dataframe(results, result_path)

    log.info(f"PubChem search completed: {len(results)} analogs found")
    return result_path
