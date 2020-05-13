########################################################################################################################

__doc__ = \
    """
These are extras.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

import logging
import os
import pymol2
import re
import requests
import sys
import unicodedata
from typing import List, Union, Optional

from rdkit import Chem
from rdkit.Chem import rdFMCS

from ._victor_base_mixin import _VictorBaseMixin


class _VictorUtilsMixin(_VictorBaseMixin):
    # =================== Logging ======================================================================================

    @classmethod
    def enable_stdout(cls, level=logging.INFO) -> None:
        """
        The ``cls.journal`` is output to the terminal.

        :param level: logging level
        :return: None
        """
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)

    @classmethod
    def enable_logfile(cls, filename='reanimation.log', level=logging.INFO) -> None:
        """
        The journal is output to a file.

        :param filename: file to write.
        :param level: logging level
        :return: None
        """
        handler = logging.FileHandler(filename)
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
        cls.journal.addHandler(handler)
        # logging.getLogger('py.warnings').addHandler(handler)

    @classmethod
    def slack_me(cls, msg: str) -> bool:
        """
        Send message to a slack webhook

        :param msg: Can be dirty and unicode-y.
        :return: did it work?
        :rtype: bool
        """
        webhook = os.environ['SLACK_WEBHOOK']
        # sanitise.
        msg = unicodedata.normalize('NFKD', msg).encode('ascii', 'ignore').decode('ascii')
        msg = re.sub('[^\w\s\-.,;?!@#()\[\]]', '', msg)
        r = requests.post(url=webhook,
                          headers={'Content-type': 'application/json'},
                          data=f"{{'text': '{msg}'}}")
        if r.status_code == 200 and r.content == b'ok':
            return True
        else:
            return False

    # =================== Other ========================================================================================

    @classmethod
    def copy_names(cls, acceptor_mol: Chem.Mol, donor_mol: Chem.Mol):
        """
        Copy names form donor to acceptor by finding MCS.
        Does it properly and uses ``PDBResidueInfo``.

        :param acceptor_mol: needs atomnames
        :param donor_mol: has atomnames
        :return:
        """
        mcs = rdFMCS.FindMCS([acceptor_mol, donor_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        pos_match = acceptor_mol.positioned_mol.GetSubstructMatch(common)
        pdb_match = donor_mol.GetSubstructMatch(common)
        for m, p in zip(pos_match, pdb_match):
            ma = acceptor_mol.GetAtomWithIdx(m)
            pa = donor_mol.GetAtomWithIdx(p)
            assert ma.GetSymbol() == pa.GetSymbol(), 'The indices do not align! ' + \
                                                     f'{ma.GetIdx()}:{ma.GetSymbol()} vs. ' + \
                                                     f'{pa.GetIdx()}:{pa.GetSymbol()}'
            ma.SetMonomerInfo(pa.GetPDBResidueInfo())

    @classmethod
    def add_constraint_to_warhead(cls, name: str, constraint: str):
        """
        Add a constraint (multiline is fine) to a warhead definition.
        This will be added and run by Egor's minimiser.

        :param name:
        :param constraint:
        :return: None
        """
        for war_def in cls.warhead_definitions:
            if war_def['name'] == name:
                war_def['constraint'] = constraint
                break
        else:
            raise ValueError(f'{name} not found in warhead_definitions.')

    @classmethod
    def closest_hit(cls, pdb_filenames: List[str],
                    target_resi: int,
                    target_chain: str,
                    target_atomname: str,
                    ligand_resn='LIG') -> str:
        """
        This classmethod helps choose which pdb based on which is closer to a given atom.

        :param pdb_filenames:
        :param target_resi:
        :param target_chain:
        :param target_atomname:
        :param ligand_resn:
        :return:
        """
        best_d = 99999
        best_hit = -1
        with pymol2.PyMOL() as pymol:
            for hit in pdb_filenames:
                pymol.cmd.load(hit)
                d = min(
                    [pymol.cmd.distance(f'chain {target_chain} and resi {target_resi} and name {target_atomname}',
                                        f'resn {ligand_resn} and name {atom.name}') for atom in
                     pymol.cmd.get_model(f'resn {ligand_resn}').atom])
                if d < best_d:
                    best_hit = hit
                    best_d = d
                pymol.cmd.delete('*')
        return best_hit

    @classmethod
    def make_covalent(cls, smiles: str, warhead_name: Optional[str] = None) -> Union[str, None]:
        """
        Convert a unreacted warhead to a reacted one in the SMILES

        :param smiles: unreacted SMILES
        :param warhead_name: name in the definitions. If unspecified it will try and guess (less preferrable)
        :return: SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        if warhead_name:
            war_defs = [wd for wd in cls.warhead_definitions if wd['name'] == warhead_name]
        else:
            war_defs = cls.warhead_definitions
        for war_def in war_defs:
            ncv = Chem.MolFromSmiles(war_def['noncovalent'])
            cv = Chem.MolFromSmiles(war_def['covalent'])
            if mol.HasSubstructMatch(ncv):
                x = Chem.ReplaceSubstructs(mol, ncv, cv, replacementConnectionPoint=0)[0]
                return Chem.MolToSmiles(x)
        else:
            return None

    @classmethod
    def make_all_warhead_combinations(cls, smiles: str, warhead_name: str, canonical=True) -> Union[dict, None]:
        """
        Convert a unreacted warhead to a reacted one in the SMILES

        :param smiles: unreacted SMILES
        :param warhead_name: name in the definitions
        :param canonical: the SMILES canonical? (makes sense...)
        :return: dictionary of SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        war_def = [wd for wd in cls.warhead_definitions if wd['name'] == warhead_name][0]
        ncv = Chem.MolFromSmiles(war_def['noncovalent'])
        if mol.HasSubstructMatch(ncv):
            combinations = {}
            for wd in cls.warhead_definitions:
                x = Chem.ReplaceSubstructs(mol, ncv, Chem.MolFromSmiles(wd['covalent']),
                                           replacementConnectionPoint=0)
                combinations[wd['name'] + '_covalent'] = Chem.MolToSmiles(x[0], canonical=canonical)
                x = Chem.ReplaceSubstructs(mol, ncv, Chem.MolFromSmiles(wd['noncovalent']),
                                           replacementConnectionPoint=0)
                combinations[wd['name'] + '_noncovalent'] = Chem.MolToSmiles(x[0], canonical=canonical)
            return combinations
        else:
            return None

    # =================== Laboratory ===================================================================================

    @classmethod
    def laboratory(cls, entries: List[dict], cores: int = 1):
        raise NotImplementedError('Not yet written.')
        pass
