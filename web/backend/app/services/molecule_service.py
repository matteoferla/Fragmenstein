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
            # Extract each ligand instance separately from the PDB.
            # A single PDB may contain multiple fragment copies (different
            # residue numbers or chains). We split them so each becomes its
            # own hit molecule for Fragmenstein.
            from fragmenstein import Victor
            pdb_text = filepath.read_text()
            instances = _find_ligand_instances(pdb_text, ligand_resn)
            base_name = filepath.stem

            for chain, resi_num in instances:
                instance_block = _extract_instance_pdb(pdb_text, ligand_resn, chain, resi_num)
                inst_name = f"{base_name}_{chain}_{resi_num}" if len(instances) > 1 else base_name
                smiles = smiles_map.get(inst_name) if smiles_map else None
                if smiles is None and smiles_map:
                    smiles = smiles_map.get(base_name)
                try:
                    mol = Victor.extract_mol(
                        name=inst_name,
                        block=instance_block,
                        smiles=smiles,
                        ligand_resn=ligand_resn,
                        proximityBonding=proximity_bonding,
                    )
                    if mol is not None:
                        mol.SetProp("_Name", inst_name)
                        mols.append(mol)
                except Exception as e:
                    log.warning(f"Failed to extract ligand '{ligand_resn}' instance {chain}/{resi_num} from {filepath.name}: {e}")
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


def _find_ligand_instances(pdb_text: str, ligand_resn: str) -> list[tuple[str, int]]:
    """Find distinct (chain, residue_number) pairs for a given residue name in a PDB."""
    seen: set[tuple[str, int]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname != ligand_resn:
            continue
        chain = line[21].strip() or "A"
        try:
            resi_num = int(line[22:26].strip())
        except ValueError:
            continue
        seen.add((chain, resi_num))
    return sorted(seen)


def _extract_instance_pdb(pdb_text: str, ligand_resn: str, chain: str, resi_num: int) -> str:
    """Extract a PDB block containing only the protein + one specific ligand instance."""
    lines = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname == ligand_resn:
                line_chain = line[21].strip() or "A"
                try:
                    line_resi = int(line[22:26].strip())
                except ValueError:
                    continue
                if line_chain != chain or line_resi != resi_num:
                    continue
        lines.append(line)
    return "\n".join(lines)


def mol_info(mol: Chem.Mol) -> dict:
    """Extract basic info from a molecule."""
    return {
        "name": mol_name(mol),
        "smiles": mol_to_smiles(mol),
        "num_atoms": mol_num_heavy_atoms(mol),
    }
