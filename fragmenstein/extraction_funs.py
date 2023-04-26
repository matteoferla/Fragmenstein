"""
These are common functions moved out of Victor and Igor to reduce duplication.
"""
import logging
import warnings

from rdkit import Chem
from typing import Tuple, Union, Callable

log = logging.getLogger()

def add_dummy_to_mol(ligand: Chem.Mol,
                     ligand_resn: str,
                     holo: Chem.Mol,
                     dummy_name: str = 'CONN') -> Chem.Mol:
    """
    Given a ligand name and a holo PDB as a Chem.Mol extract the dummy atom and add it to the ligand Chem.Mol.
    Used by Victor and Igor in different points.

    :param ligand:
    :param ligand_resn:
    :param holo:
    :return:
    """
    attachment, attachee = find_attachment(holo, ligand_resn)
    if attachment is not None:  # covalent
        mod = Chem.RWMol(ligand)
        attachment.SetAtomicNum(0)  # dummy atom.
        attachment.GetPDBResidueInfo().SetName(dummy_name)
        pos = holo.GetConformer().GetAtomPosition(attachment.GetIdx())
        ni = mod.AddAtom(attachment)
        mod.GetConformer().SetAtomPosition(ni, pos)
        attachee_name = attachee.GetPDBResidueInfo().GetName()
        for atom in mod.GetAtoms():
            if atom.GetPDBResidueInfo().GetName() == attachee_name:
                ai = atom.GetIdx()
                mod.AddBond(ai, ni, Chem.BondType.SINGLE)
                break
        ligand = mod.GetMol()
        ligand.UpdatePropertyCache(strict=False)
    return ligand


def find_attachment(pdb: Chem.Mol, ligand_resn: str) -> Tuple[Union[Chem.Atom, None], Union[Chem.Atom, None]]:
    """
    Finds the two atoms in a crosslink bond without looking at LINK record

    :param pdb: a rdkit Chem object
    :param ligand_resn: 3 letter code
    :return: tuple of non-ligand atom and ligand atom
    """
    for atom in pdb.GetAtoms():
        if atom.GetPDBResidueInfo().GetResidueName() == ligand_resn:
            for neigh in atom.GetNeighbors():
                if neigh.GetPDBResidueInfo().GetResidueName() != ligand_resn:
                    attachment = neigh
                    attachee = atom
                    return attachment, attachee
    attachment = None
    attachee = None
    return attachment, attachee


def get_atomname(atom: Chem.Atom) -> str:
    return atom.GetPDBResidueInfo().GetName()

def get_equivalent(template_atom: Chem.Atom, query_mol: Chem.Mol) -> Chem.Atom:
    """
    Given an atom find the same PDB named one in the query mol.
    """
    get_by_atomname: Callable[[Chem.Mol, str], Chem.Atom] = \
        lambda mol, name: [atom for atom in mol.GetAtoms() if get_atomname(atom)[:4] == name[:4]][0]
    return get_by_atomname(query_mol, get_atomname(template_atom))

def combine_for_bondorder(template: Chem.Mol, target: Chem.Mol) -> Chem.Mol:
    target2template = {}
    for template_ai, template_atom in enumerate(template.GetAtoms()):
        target_atom: Chem.Atom = get_equivalent(template_atom, target)
        target_ai: int = target_atom.GetIdx()
        assert target_ai not in target2template.keys(), f'The atom name {get_atomname(target_atom)} is duplicated '+\
                                                   f'in the template {get_atomname(template_atom)} & '+\
                                                   f'{get_atomname(template.GetAtomWithIdx(target2template[target_ai]))}'
        target2template[target_ai] = template_ai
    new = Chem.Mol(template)
    # fix coordinates
    conf: Chem.Conformer = new.GetConformer()
    donor_conf: Chem.Conformer = target.GetConformer()
    for target_ai, template_ai in target2template.items():
        conf.SetAtomPosition(template_ai, donor_conf.GetAtomPosition(target_ai))
    # fix props
    for prop, value in target.GetPropsAsDict().items():
        if not new.HasProp(prop):
            new.SetProp(prop, value)
    for target_ai, new_ai in target2template.items():
        target_atom = target.GetAtomWithIdx(target_ai)
        new_atom = new.GetAtomWithIdx(new_ai)
        for prop, value in target_atom.GetPropsAsDict().items():
            if new_atom.HasProp(prop):
                continue
            new_atom.SetProp(prop, value)
    return new


def copy_bonds_by_atomnames(template: Chem.Mol, target: Chem.Mol) -> bool:
    """
    Fixes bonds _inplace_ and sanity checks the mol is still the same.

    ``AllChem.AssignBondOrdersFromTemplate(template, target)``
    may not respect the atoms. Say cyclohexane -> cyclohexene alone can be mapped 12 ways by MCS,
    but only one by atomname.

    Plus ``PDBResidueInfo`` is lost by ``AllChem.AssignBondOrdersFromTemplate``
    """
    warnings.warn('This is not used anymore, but kept for reference', DeprecationWarning)
    successful = True
    for template_bond in template.GetBonds():  # type: Chem.Bond
        begin: Chem.Atom = template_bond.GetBeginAtom()
        end: Chem.Atom = template_bond.GetEndAtom()
        query_bond = target.GetBondBetweenAtoms(get_equivalent(begin, target).GetIdx(),
                                                get_equivalent(end, target).GetIdx())
        if query_bond is None:
            # This is due to the molecule being pulled apart at a ring closure.
            # A common issue with RNA minisation too.
            print(f'{get_atomname(begin)} and {get_atomname(end)} are' + \
                             ' bonded in the params but not the minimised: this is a Rosetta glitch')
            successful = False
            continue
        query_bond.SetBondType(template_bond.GetBondType())
        # heme will have crash upstream though:
        query_bond.SetBondDir(template_bond.GetBondDir())
    for template_atom in template.GetAtoms():
        target_atom: Chem.Atom = get_equivalent(template_atom, target)
        #print(get_atomname(template_atom), get_atomname(target_atom))
        target_atom.SetIsAromatic(template_atom.GetIsAromatic())
        target_atom.SetIsotope(template_atom.GetIsotope())
        target_atom.SetFormalCharge(template_atom.GetFormalCharge())
        target_atom.SetChiralTag(template_atom.GetChiralTag())
        target_atom.SetNumRadicalElectrons(template_atom.GetNumRadicalElectrons())
        target_atom.SetNumExplicitHs(template_atom.GetNumExplicitHs())
        target_atom.SetHybridization(template_atom.GetHybridization())
        target_atom.SetAtomMapNum(template_atom.GetAtomMapNum())
    # problem here:
    target.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(target)
    return successful
