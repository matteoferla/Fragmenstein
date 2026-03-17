"""Single-molecule Victor operations — run one combine or place."""

import math
import logging

from rdkit import Chem

from . import file_manager
from .combine_service import VICTOR_TYPES, load_session_hits, load_session_template

log = logging.getLogger(__name__)


def _safe_float(v):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return None
    return v


def run_single_combine(
    session_id: str,
    hit_names: list[str],
    victor_type: str = "Wictor",
    warhead_harmonisation: str = "first",
    joining_cutoff: float = 5.0,
    covalent_resi: str | None = None,
) -> dict:
    """Run Victor.combine() on a single set of hits."""
    hits = load_session_hits(session_id, hit_names)
    if len(hits) < 2:
        raise ValueError("Need at least 2 hits")

    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)

    VictorClass = VICTOR_TYPES.get(victor_type)
    if VictorClass is None:
        from fragmenstein import Wictor
        VictorClass = Wictor

    VictorClass.work_path = str(work_dir)
    VictorClass.monster_joining_cutoff = joining_cutoff
    VictorClass.error_to_catch = Exception

    victor = VictorClass(
        hits=hits,
        pdb_block=pdbblock,
        ligand_resn="LIG",
        covalent_resi=covalent_resi,
    )
    victor.monster_throw_on_discard = True
    victor.combine(warhead_harmonisation=warhead_harmonisation)

    summary = victor.summarize()
    result = {
        "name": summary.get("name", ""),
        "smiles": summary.get("smiles"),
        "error": summary.get("error"),
        "mode": summary.get("mode"),
        "ddG": _safe_float(summary.get("∆∆G")),
        "dG_bound": _safe_float(summary.get("∆G_bound")),
        "dG_unbound": _safe_float(summary.get("∆G_unbound")),
        "comRMSD": _safe_float(summary.get("comRMSD")),
        "N_constrained_atoms": summary.get("N_constrained_atoms"),
        "N_unconstrained_atoms": summary.get("N_unconstrained_atoms"),
        "runtime": _safe_float(summary.get("runtime")),
        "mol_block": Chem.MolToMolBlock(victor.minimized_mol) if victor.minimized_mol else None,
        "unmin_mol_block": Chem.MolToMolBlock(victor.monster.positioned_mol) if victor.monster.positioned_mol else None,
    }
    return result


def run_single_place(
    session_id: str,
    smiles: str,
    hit_names: list[str],
    victor_type: str = "Wictor",
    merging_mode: str = "expansion",
    covalent_resi: str | None = None,
) -> dict:
    """Run Victor.place() on a single molecule."""
    hits = load_session_hits(session_id, hit_names)
    if not hits:
        raise ValueError("No hits loaded")

    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)

    VictorClass = VICTOR_TYPES.get(victor_type)
    if VictorClass is None:
        from fragmenstein import Wictor
        VictorClass = Wictor

    VictorClass.work_path = str(work_dir)
    VictorClass.error_to_catch = Exception

    victor = VictorClass(
        hits=hits,
        pdb_block=pdbblock,
        ligand_resn="LIG",
        covalent_resi=covalent_resi,
    )
    victor.place(smiles, merging_mode=merging_mode)

    summary = victor.summarize()
    result = {
        "name": summary.get("name", "placed"),
        "smiles": summary.get("smiles"),
        "error": summary.get("error"),
        "mode": summary.get("mode"),
        "ddG": _safe_float(summary.get("∆∆G")),
        "dG_bound": _safe_float(summary.get("∆G_bound")),
        "dG_unbound": _safe_float(summary.get("∆G_unbound")),
        "comRMSD": _safe_float(summary.get("comRMSD")),
        "N_constrained_atoms": summary.get("N_constrained_atoms"),
        "N_unconstrained_atoms": summary.get("N_unconstrained_atoms"),
        "runtime": _safe_float(summary.get("runtime")),
        "mol_block": Chem.MolToMolBlock(victor.minimized_mol) if victor.minimized_mol else None,
        "unmin_mol_block": Chem.MolToMolBlock(victor.monster.positioned_mol) if victor.monster.positioned_mol else None,
    }
    return result
