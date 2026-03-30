"""Monster standalone operations — chemistry-only, no protein context."""

import logging

from rdkit import Chem

from fragmenstein import Monster

from .combine_service import load_session_hits

log = logging.getLogger(__name__)


def run_monster_combine(
    session_id: str,
    joining_cutoff: float = 5.0,
    hit_names: list[str] | None = None,
) -> dict:
    """Run Monster.combine() without protein. Returns result dict."""
    hits = load_session_hits(session_id, hit_names)
    if len(hits) < 2:
        raise ValueError("Need at least 2 hits for combine")

    monster = Monster(hits=hits)
    monster.joining_cutoff = joining_cutoff
    monster.combine()

    mol = monster.positioned_mol
    name = "-".join(m.GetProp("_Name") for m in hits if m.HasProp("_Name"))

    return {
        "name": name,
        "smiles": Chem.MolToSmiles(mol) if mol else None,
        "num_atoms": mol.GetNumHeavyAtoms() if mol else None,
        "mol_block": Chem.MolToMolBlock(mol) if mol else None,
    }


def run_monster_place(
    session_id: str,
    smiles: str,
    hit_names: list[str] | None = None,
) -> dict:
    """Run Monster.place() without protein. Returns result dict."""
    hits = load_session_hits(session_id, hit_names)
    if not hits:
        raise ValueError("No hits loaded")

    monster = Monster(hits=hits)
    monster.place(smiles)

    mol = monster.positioned_mol

    return {
        "name": "placed",
        "smiles": Chem.MolToSmiles(mol) if mol else None,
        "num_atoms": mol.GetNumHeavyAtoms() if mol else None,
        "mol_block": Chem.MolToMolBlock(mol) if mol else None,
    }
