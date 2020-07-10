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
import re
import requests
import sys, json
import unicodedata
from typing import List, Union, Optional, Dict, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem

from ._victor_base_mixin import _VictorBaseMixin
from ..m_rmsd import mRSMD
from ..core import Fragmenstein
from rdkit_to_params import Params
from ..igor import Igor

try:
    import pymol2
except ImportError:
    pymol2 = None


class _VictorUtilsMixin(_VictorBaseMixin):

    def dock(self) -> Chem.Mol:
        docked = self.igor.dock()
        self.docked_pose = docked
        docked.dump_pdb(f'{self.work_path}/{self.long_name}/{self.long_name}.holo_docked.pdb')
        ligand = self.igor.mol_from_pose(docked)
        template = AllChem.DeleteSubstructs(self.params.mol, Chem.MolFromSmiles('*'))
        lig_chem = AllChem.AssignBondOrdersFromTemplate(template, ligand)
        Chem.MolToMolFile(lig_chem, f'{self.work_path}/{self.long_name}/{self.long_name}.docked.mol')
        return lig_chem
        # print(pyrosetta.get_fa_scorefxn()(docked) - v.energy_score['unbound_ref2015']['total_score'])

    def summarise(self):
        return {'name': self.long_name,
                'smiles': self.smiles,
                'mode': self.fragmenstein_merging_mode,
                '∆∆G': self.energy_score['ligand_ref2015']['total_score'] - \
                       self.energy_score['unbound_ref2015']['total_score'],
                '∆G_bound': self.energy_score['ligand_ref2015']['total_score'],
                '∆G_unbound': self.energy_score['unbound_ref2015']['total_score'],
                'comRMSD': self.mrmsd.mrmsd,
                'N_constrained_atoms': self.constrained_atoms,
                'N_unconstrained_atoms': self.unconstrained_heavy_atoms,
                'runtime': self.tock - self.tick,
                'regarded': [h.GetProp('_Name') for h in self.hits if
                             h.GetProp('_Name') not in self.fragmenstein.unmatched],
                'disregarded': self.fragmenstein.unmatched
                }

    # =================== Logging ======================================================================================

    @classmethod
    def enable_stdout(cls, level=logging.INFO) -> None:
        """
        The ``cls.journal`` is output to the terminal.

        :param level: logging level
        :return: None
        """
        cls.journal.handlers = [h for h in cls.journal.handlers if h.name != 'stdout']
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(level)
        handler.set_name('stdout')
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
        cls.journal.handlers = [h for h in cls.journal.handlers if h.name != 'logfile']
        handler = logging.FileHandler(filename)
        handler.setLevel(level)
        handler.set_name('logfile')
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
        This will be added and run by Igor's minimiser.

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
    def distance_hits(cls, pdb_filenames: List[str],
                      target_resi: int,
                      target_chain: str,
                      target_atomname: str,
                      ligand_resn='LIG') -> List[float]:
        """
        See closest hit for info.

        :param pdb_filenames:
        :param target_resi:
        :param target_chain:
        :param target_atomname:
        :param ligand_resn:
        :return:
        """
        distances = []
        with pymol2.PyMOL() as pymol:
            for hit in pdb_filenames:
                pymol.cmd.load(hit)
                distances.append(min(
                    [pymol.cmd.distance(f'chain {target_chain} and resi {target_resi} and name {target_atomname}',
                                        f'resn {ligand_resn} and name {atom.name}') for atom in
                     pymol.cmd.get_model(f'resn {ligand_resn}').atom]))
                pymol.cmd.delete('*')
        return distances

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
        for hit, d in zip(pdb_filenames,
                          cls.distance_hits(pdb_filenames, target_resi, target_chain, target_atomname, ligand_resn)):
            if d < best_d:
                best_hit = hit
                best_d = d
        return best_hit

    @classmethod
    def make_covalent(cls, smiles: str,
                      warhead_name: Optional[str] = None) -> Union[str, None]:
        """
        Convert a unreacted warhead to a reacted one in the SMILES

        :param smiles: unreacted SMILES
        :param warhead_name: name in the definitions. If unspecified it will try and guess (less preferrable)
        :return: SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        if warhead_name:
            war_defs = cls._get_warhead_definitions(warhead_name)
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
    def get_warhead_definition(cls, warhead_name: str):
        return cls._get_warhead_definitions(warhead_name)[0]

    @classmethod
    def _get_warhead_definitions(cls, warhead_name: str):
        """
        It is unlikely that alternative definitions are present. hence why hidden method.

        :param warhead_name:
        :return:
        """
        options = [wd for wd in cls.warhead_definitions if wd['name'] == warhead_name.lower()]
        if len(options) == 0:
            raise ValueError(f'{warhead_name} is not valid.')
        else:
            return options

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
        war_def = cls.get_warhead_definition(warhead_name)
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

    # =================== pre-encounter ================================================================================

    # @classmethod

    # =================== save  ========================================================================================

    def make_pse(self, filename: str = 'combo.pse'):
        """
        Save a pse in the relevant folder.

        :param filename:
        :return:
        """
        assert '.pse' in filename, f'{filename} not .pse file'
        with pymol2.PyMOL() as pymol:
            for hit in self.hits:
                hit_name = hit.GetProp('_Name')
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit), hit_name)
                if hit_name in self.fragmenstein.unmatched:
                    pymol.cmd.color('black', f'element C and {hit_name}')
                else:
                    pymol.cmd.color('white', f'element C and {hit_name}')
            if self.fragmenstein.positioned_mol is not None:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.fragmenstein.positioned_mol), 'placed')
                pymol.cmd.color('magenta', f'element C and placed')
            if self.minimised_mol is not None:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.minimised_mol), 'minimised')
                pymol.cmd.color('green', f'element C and minimised')
            if self.minimised_pdbblock is not None:
                pymol.cmd.read_pdbstr(self.minimised_pdbblock, 'min_protein')
                pymol.cmd.color('gray50', f'element C and min_protein')
                pymol.cmd.hide('sticks', 'min_protein')
            if self.unminimised_pdbblock is not None:
                pymol.cmd.read_pdbstr(self.unminimised_pdbblock, 'unmin_protein')
                pymol.cmd.color('gray20', f'element C and unmin_protein')
                pymol.cmd.hide('sticks', 'unmin_protein')
                pymol.cmd.disable('unmin_protein')
            pymol.cmd.zoom('byres (placed expand 4)')
            pymol.cmd.show('line', 'byres (placed around 4)')
            pymol.cmd.save(os.path.join(self.work_path, self.long_name, filename))

    # =================== extract_mols =================================================================================

    @classmethod
    def find_attachment(self, pdb: Chem.Mol, ligand_resn: str) -> Tuple[Union[Chem.Atom, None], Union[Chem.Atom, None]]:
        for atom in pdb.GetAtoms():
            if atom.GetPDBResidueInfo().GetResidueName() == ligand_resn:
                for neigh in atom.GetNeighbors():
                    if neigh.GetPDBResidueInfo().GetResidueName() != ligand_resn:
                        attachment = neigh
                        attachee = atom
                        return (attachment, attachee)
        else:
            attachment = None
            attachee = None
            return (attachment, attachee)

    @classmethod
    def extract_mols(cls,
                     folder: str,
                     smilesdex: Dict[str, str],
                     ligand_resn: str = 'LIG',
                     regex_name: Optional[str]= None) -> Dict[str, Chem.Mol]:
        """
         A key requirement for Fragmenstein is a separate mol file for the inspiration hits.
        This is however often a pdb. This converts.
        `igor.mol_from_pose()` is similar but works on a pose. `_fix_minimised()` calls mol_from_pose.

        :param folder: folder with pdbs
        :return:
        """
        mols = {}
        for file in os.listdir(folder):
            if '.pdb' not in file:
                continue
            else:
                fullfile = os.path.join(folder, file)
                if regex_name is None:
                    name = os.path.splitext(file)[0]
                elif re.search(regex_name, file) is None:
                    continue
                else:
                    name = re.search(regex_name, file).group(1)
                if name in smilesdex:
                    smiles=smilesdex[name]
                else:
                    cls.journal.warning(f'{name} could not be matched to a smiles.')
                    smiles = None
                try:
                    mol = cls.extract_mol(name=name,
                                          filepath=fullfile,
                                          smiles=smiles,
                                          ligand_resn=ligand_resn)
                    if mol is not None:
                        mols[name] = mol
                except Exception as error:
                    cls.journal.error(f'{error.__class__.__name__} for {name} - {error}')
        return mols

    @classmethod
    def extract_mol(cls,
                     name: str,
                     filepath: str,
                     smiles: Optional[str] = None,
                     ligand_resn: str = 'LIG') -> Chem.Mol:
        holo = Chem.MolFromPDBFile(filepath, proximityBonding=False, removeHs=False)
        mol = Chem.SplitMolByPDBResidues(holo, whiteList=[ligand_resn])[ligand_resn]
        attachment, attachee = cls.find_attachment(holo, ligand_resn)
        if attachment is not None:  # covalent
            mol = Chem.SplitMolByPDBResidues(holo, whiteList=[ligand_resn])[ligand_resn]
            mod = Chem.RWMol(mol)
            attachment.SetAtomicNum(0)  # dummy atom.
            attachment.GetPDBResidueInfo().SetName('CONN')
            pos = holo.GetConformer().GetAtomPosition(attachment.GetIdx())
            ni = mod.AddAtom(attachment)
            mod.GetConformer().SetAtomPosition(ni, pos)
            attachee_name = attachee.GetPDBResidueInfo().GetName()
            for atom in mod.GetAtoms():
                if atom.GetPDBResidueInfo().GetName() == attachee_name:
                    ai = atom.GetIdx()
                    mod.AddBond(ai, ni, Chem.BondType.SINGLE)
                    break
            mol = mod.GetMol()
        if smiles is not None:
            if '*' in Chem.MolToSmiles(mol) and '*' not in smiles:
                smiles = cls.make_covalent(smiles)
                cls.journal.info(f'{name} is covalent but a non covalent SMILES was passed.')
            else:
                pass
            try:
                template = Chem.MolFromSmiles(smiles)
                # template = AllChem.DeleteSubstructs(template, Chem.MolFromSmiles('*'))
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except ValueError as error:
                cls.journal.warning(f'{name} failed at bonding ({type(error)}: {error}).')
        mol.SetProp('_Name', name)
        return mol








    # =================== From Files ===================================================================================

    @classmethod
    def from_files(cls, folder: str) -> _VictorBaseMixin:
        """
        This creates an instance form the output files. Likely to be unstable.
        Assumes the checkpoints were not altered.
        And is basically for analysis only.

        :param folder: path
        :return:
        """
        cls.journal.warning('`from_files`: You really should not use this.')
        if os.path.exists(folder):
            pass # folder is fine
        elif not os.path.exists(folder) and os.path.exists(os.path.join(cls.work_path, folder)):
            folder = os.path.join(cls.work_path, folder)
        else:
            raise FileNotFoundError(f'Folder {folder} does not exist.')
        self = cls.__new__(cls)
        self.tick = float('nan')
        self.tock = float('nan')
        self.ligand_resn = ''
        self.ligand_resi = ''
        self.covalent_resn = ''
        self.covalent_resi = ''
        self.hits = []
        self.long_name = os.path.split(folder)[1]
        paramsfiles = os.path.join(folder, f'{self.long_name}.params')
        paramstemp = os.path.join(folder, f'{self.long_name}.params_template.mol')
        if os.path.exists(paramsfiles):
            self.params = Params().load(paramsfiles)
            self.unbound_pose = self.params.test()
            if os.path.exists(paramstemp):
                self.params.mol = Chem.MolFromMolFile(paramstemp, removeHs=False)
        else:
            self.params = None
        posmol = os.path.join(folder, f'{self.long_name}.positioned.mol')
        if os.path.exists(posmol):
            self.mol = Chem.MolFromMolFile(posmol, sanitize=False, removeHs=False)
        else:
            self.journal.info(f'{self.long_name} - no positioned mol')
            self.mol = None
        fragjson = os.path.join(folder, f'{self.long_name}.fragmenstein.json')
        if os.path.exists(fragjson):
            fd = json.load(open(fragjson))
            self.smiles = fd['smiles']
            self.is_covalent = True if '*' in self.smiles else False
            self.fragmenstein = Fragmenstein(mol=self.mol,
                                             hits=self.hits,
                                             attachment=None,
                                             merging_mode='off',
                                             average_position=self.fragmenstein_average_position
                                             )
            self.fragmenstein.positioned_mol = self.mol
            self.fragmenstein.positioned_mol.SetProp('_Origins', json.dumps(fd['origin']))

        else:
            self.is_covalent = None
            self.smiles = ''
            self.fragmenstein = None
            self.journal.info(f'{self.long_name} - no fragmenstein json')
            self.N_constrained_atoms = float('nan')

        #
        self.apo_pdbblock = None
        #
        self.atomnames = None
        #
        self.extra_constraint = ''
        self.pose_fx = None
        # these are calculated
        self.constraint = None
        self.unminimised_pdbblock = None
        self.igor = None
        self.minimised_pdbblock = None
        # buffers etc.
        self._warned = []
        minjson = os.path.join(folder, f'{self.long_name}.minimised.json')
        self.mrmsd = mRSMD.mock()
        if os.path.exists(minjson):
            md = json.load(open(minjson))
            self.energy_score = md["Energy"]
            self.mrmsd.mrmsd = md["mRMSD"]
            self.mrmsd.rmsds = md["RMSDs"]
            self.igor = Igor.from_pdbfile(
                pdbfile=os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb'),
                params_file=os.path.join(self.work_path, self.long_name, self.long_name + '.params'),
                constraint_file=os.path.join(self.work_path, self.long_name, self.long_name + '.con'))
        else:
            self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                                 'unbound_ref2015': {'total_score': float('nan')}}

            self.journal.info(f'{self.long_name} - no min json')
        minmol = os.path.join(folder, f'{self.long_name}.minimised.mol')
        if os.path.exists(minmol):
            self.minimised_mol = Chem.MolFromMolFile(minmol, sanitize=False, removeHs=False)
        else:
            self.minimised_mol = None
        return self

    # =================== Laboratory ===================================================================================

    @classmethod
    def laboratory(cls, entries: List[dict], cores: int = 1):
        raise NotImplementedError('Not yet written.')
        pass
