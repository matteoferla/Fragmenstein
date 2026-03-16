"""RDKit molecule parsing — SDF/MOL/PDB to Mol objects."""

import logging
from pathlib import Path

from rdkit import Chem

log = logging.getLogger(__name__)


def parse_smi_file(filepath: Path) -> dict[str, str]:
    """Parse a .smi file returning a name→SMILES mapping.

    Format: SMILES<tab>NAME per line.
    """
    mapping = {}
    text = filepath.read_text()
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            smiles, name = parts[0].strip(), parts[-1].strip()
            if name and smiles:
                mapping[name] = smiles
    return mapping


def parse_mol_file(
    filepath: Path,
    ligand_resn: str | None = None,
    proximity_bonding: bool = True,
    smiles_map: dict[str, str] | None = None,
) -> list[Chem.Mol]:
    """Parse a molecular file and return a list of RDKit Mol objects.

    Supports .sdf, .mol, .mol2, .pdb formats.
    For SDF files, may return multiple molecules.
    For PDB files: if ligand_resn is provided, extracts only the ligand
    residue using Victor.extract_mol (standard for crystal structure hits).
    If not provided, reads the whole PDB as a single molecule.
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
        if ligand_resn:
            # Extract ligand from crystal structure PDB using Victor
            from fragmenstein import Victor
            name = filepath.stem
            # Look up SMILES for bond order correction
            smiles = smiles_map.get(name) if smiles_map else None
            try:
                mol = Victor.extract_mol(
                    name=name,
                    filepath=str(filepath),
                    smiles=smiles,
                    ligand_resn=ligand_resn,
                    proximityBonding=proximity_bonding,
                )
                if mol is not None:
                    mol.SetProp("_Name", name)
                    mols.append(mol)
            except Exception as e:
                log.warning(f"Failed to extract ligand '{ligand_resn}' from {filepath.name}: {e}")
        else:
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
