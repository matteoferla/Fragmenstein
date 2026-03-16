"""RDKit molecule parsing — SDF/MOL/PDB to Mol objects."""

import logging
from pathlib import Path

from rdkit import Chem

log = logging.getLogger(__name__)


def parse_mol_file(filepath: Path) -> list[Chem.Mol]:
    """Parse a molecular file and return a list of RDKit Mol objects.

    Supports .sdf, .mol, .mol2, .pdb formats.
    For SDF files, may return multiple molecules.
    """
    suffix = filepath.suffix.lower()
    mols = []

    if suffix == ".sdf":
        supplier = Chem.SDMolSupplier(str(filepath), removeHs=False)
        for mol in supplier:
            if mol is not None:
                mols.append(mol)
    elif suffix in (".mol", ".mol2"):
        mol = Chem.MolFromMolFile(str(filepath), removeHs=False)
        if mol is not None:
            mols.append(mol)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(filepath), removeHs=False)
        if mol is not None:
            mols.append(mol)
    else:
        raise ValueError(f"Unsupported file format: {suffix}")

    return mols


def mol_to_smiles(mol: Chem.Mol) -> str | None:
    """Convert an RDKit Mol to SMILES string."""
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def mol_to_molblock(mol: Chem.Mol) -> str:
    """Convert an RDKit Mol to MolBlock string (for 3D viewer)."""
    return Chem.MolToMolBlock(mol)


def mol_name(mol: Chem.Mol) -> str:
    """Get the name of a molecule from its properties."""
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    return name or "unnamed"


def mol_num_heavy_atoms(mol: Chem.Mol) -> int:
    """Count heavy (non-hydrogen) atoms."""
    return mol.GetNumHeavyAtoms()


def mol_info(mol: Chem.Mol) -> dict:
    """Extract basic info from a molecule."""
    return {
        "name": mol_name(mol),
        "smiles": mol_to_smiles(mol),
        "num_atoms": mol_num_heavy_atoms(mol),
    }
