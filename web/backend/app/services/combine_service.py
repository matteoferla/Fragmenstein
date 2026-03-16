"""Wraps Laboratory.combine() for the web API."""

import logging
from pathlib import Path

from rdkit import Chem

from fragmenstein import Laboratory, Victor, Wictor

from . import file_manager
from .molecule_service import parse_mol_file
from .result_serializer import save_dataframe

log = logging.getLogger(__name__)

VICTOR_TYPES = {
    "Victor": Victor,
    "Wictor": Wictor,
}

# Import optional variants
try:
    from fragmenstein import Quicktor
    VICTOR_TYPES["Quicktor"] = Quicktor
except ImportError:
    pass

try:
    from fragmenstein import OpenVictor
    VICTOR_TYPES["OpenVictor"] = OpenVictor
except ImportError:
    pass


def load_session_hits(session_id: str, hit_names: list[str] | None = None) -> list[Chem.Mol]:
    """Load all hit molecules for a session from the database."""
    from ..database import get_db

    with get_db() as db:
        if hit_names:
            placeholders = ",".join("?" * len(hit_names))
            rows = db.execute(
                f"SELECT name, mol_block FROM hits WHERE session_id = ? AND name IN ({placeholders})",
                [session_id] + hit_names,
            ).fetchall()
        else:
            rows = db.execute(
                "SELECT name, mol_block FROM hits WHERE session_id = ?",
                (session_id,),
            ).fetchall()

    mols = []
    for row in rows:
        mol = Chem.MolFromMolBlock(row["mol_block"], removeHs=False)
        if mol is not None:
            mol.SetProp("_Name", row["name"])
            mols.append(mol)
    return mols


def load_session_template(session_id: str) -> str:
    """Load the template PDB block for a session."""
    from ..database import get_db

    with get_db() as db:
        row = db.execute(
            "SELECT template_pdb FROM sessions WHERE id = ?",
            (session_id,),
        ).fetchone()

    if row is None or row["template_pdb"] is None:
        raise ValueError(f"No template uploaded for session {session_id}")
    return row["template_pdb"]


def run_combine(
    session_id: str,
    victor_type: str = "Wictor",
    n_cores: int = -1,
    timeout: int = 240,
    combination_size: int = 2,
    permute: bool = True,
    joining_cutoff: float = 5.0,
    quick_reanimation: bool = False,
    covalent_resi: str | None = None,
    hit_names: list[str] | None = None,
) -> Path:
    """Run Laboratory.combine() and save results. Returns path to results file."""
    hits = load_session_hits(session_id, hit_names)
    if not hits:
        raise ValueError("No hits loaded for session")

    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)

    # Resolve Victor type
    VictorClass = VICTOR_TYPES.get(victor_type, Wictor)

    # Configure Victor class attributes
    Laboratory.Victor = VictorClass
    Laboratory.Victor.work_path = str(work_dir)
    Laboratory.Victor.monster_throw_on_discard = True
    Laboratory.Victor.monster_joining_cutoff = joining_cutoff
    Laboratory.Victor.quick_reanimation = quick_reanimation
    Laboratory.Victor.error_to_catch = Exception

    lab = Laboratory(pdbblock=pdbblock, covalent_resi=covalent_resi)
    df = lab.combine(
        hits,
        n_cores=n_cores,
        timeout=timeout,
        combination_size=combination_size,
        permute=permute,
    )

    # Save results
    result_path = file_manager.results_dir(session_id, "combine") / "results.pkl"
    save_dataframe(df, result_path)

    log.info(f"Combine completed for session {session_id}: {len(df)} results")
    return result_path
