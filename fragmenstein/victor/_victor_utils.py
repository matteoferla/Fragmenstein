from __future__ import annotations
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
import numpy as np
from ..extraction_funs import add_dummy_to_mol
from typing import List, Union, Optional, Dict, Tuple, TYPE_CHECKING

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, EnumerateStereoisomers
from rdkit.Chem.MolStandardize import rdMolStandardize

from ._victor_common import _VictorCommon
from ._victor_show import _VictorShow
from ..m_rmsd import mRMSD
from ..monster import Monster
from rdkit_to_params import Params
from ..igor import Igor

if TYPE_CHECKING or 'sphinx' in sys.modules:
    import pandas as pd


class _VictorUtils(_VictorShow):  # _VictorCommon -> _VictorShow

    def dock(self) -> Chem.Mol:
        """
        The docking is done by ``igor.dock()``. This basically does that, extacts ligand, saves etc.

        :return:
        """
        docked = self.igor.dock()
        self.docked_pose = docked
        docked.dump_pdb(f'{self.work_path}/{self.long_name}/{self.long_name}.holo_docked.pdb')
        ligand = self.igor.mol_from_pose(docked)
        template = AllChem.DeleteSubstructs(self.params.mol, Chem.MolFromSmiles('*'))
        lig_chem = AllChem.AssignBondOrdersFromTemplate(template, ligand)
        lig_chem.SetProp('_Name', 'docked')
        Chem.MolToMolFile(lig_chem, f'{self.work_path}/{self.long_name}/{self.long_name}.docked.mol')
        return lig_chem
        # print(pyrosetta.get_fa_scorefxn()(docked) - v.energy_score['unbound_ref2015']['total_score'])

    def summarize(self):
        if self.error_msg:
            if self.monster is None:
                N_constrained_atoms = float('nan')
                N_unconstrained_atoms = float('nan')
            elif self.monster.positioned_mol is None:
                N_constrained_atoms = float('nan')
                N_unconstrained_atoms = float('nan')
            else:
                N_constrained_atoms = self.constrained_atoms
                N_unconstrained_atoms = self.unconstrained_heavy_atoms
            return {'name': self.long_name,
                    'smiles': self.smiles,
                    'error': self.error_msg,
                    'mode': self.merging_mode,
                    '∆∆G': float('nan'),
                    '∆G_bound': float('nan'),
                    '∆G_unbound': float('nan'),
                    'comRMSD': float('nan'),
                    'N_constrained_atoms': N_constrained_atoms,
                    'N_unconstrained_atoms': N_unconstrained_atoms,
                    'runtime': self.tock - self.tick,
                    'regarded': self.monster.matched,
                    'disregarded': self.monster.unmatched
                    }
        else:
            return {'name': self.long_name,
                    'smiles': self.smiles,
                    'error': self.error_msg,
                    'mode': self.merging_mode,
                    '∆∆G': self.energy_score['ligand_ref2015']['total_score'] - \
                           self.energy_score['unbound_ref2015']['total_score'],
                    '∆G_bound': self.energy_score['ligand_ref2015']['total_score'],
                    '∆G_unbound': self.energy_score['unbound_ref2015']['total_score'],
                    'comRMSD': self.mrmsd.mrmsd,
                    'N_constrained_atoms': self.constrained_atoms,
                    'N_unconstrained_atoms': self.unconstrained_heavy_atoms,
                    'runtime': self.tock - self.tick,
                    'regarded': self.monster.matched,
                    'disregarded': self.monster.unmatched
                    }

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
        import pymol2
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

    # =================== extract_mols =================================================================================

    @classmethod
    def find_closest_to_ligand(cls, pdb: Chem.Mol, ligand_resn: str) -> Tuple[Chem.Atom, Chem.Atom]:
        """
        Find the closest atom to the ligand
        Warning requires the protein to be loaded as an rdkit.Chem.Mol

        :param pdb: a rdkit Chem.Mol object
        :param ligand_resn: 3 letter code
        :return: tuple of non-ligand atom and ligand atom
        """
        ligand = [atom.GetIdx() for atom in pdb.GetAtoms() if atom.GetPDBResidueInfo().GetResidueName() == ligand_resn]
        dm = Chem.Get3DDistanceMatrix(pdb)
        mini = np.take(dm, ligand, 0)
        mini[mini == 0] = np.nan
        mini[:, ligand] = np.nan
        a, b = np.where(mini == np.nanmin(mini))
        lig_atom = pdb.GetAtomWithIdx(ligand[int(a[0])])
        nonlig_atom = pdb.GetAtomWithIdx(int(b[0]))
        return (nonlig_atom, lig_atom)


    @classmethod
    def extract_mols(cls,
                     folder: str,
                     smilesdex: Dict[str, str],
                     ligand_resn: str = 'LIG',
                     regex_name: Optional[str]= None,
                     proximityBonding: bool = False,
                     throw_on_error:bool=False) -> Dict[str, Chem.Mol]:
        """
         A key requirement for Monster is a separate mol file for the inspiration hits.
        This is however often a pdb. This converts.
        `igor.mol_from_pose()` is similar but works on a pose. `_fix_minimized()` calls ``mol_from_pose``
        and ``copy_bonds_by_atomnames`` which does not destroy pdbinfo.
        The latter is glitchy. Use ``combine_for_bondorder``.

        See ``extract_mol`` for single.

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
                elif throw_on_error:
                    raise ValueError(f'{name} could not be matched to a smiles.')
                else:
                    cls.journal.warning(f'{name} could not be matched to a smiles.')
                    smiles = None
                try:
                    mol = cls.extract_mol(name=name,
                                          filepath=fullfile,
                                          smiles=smiles,
                                          ligand_resn=ligand_resn,
                                          proximityBonding=proximityBonding)
                    if mol is not None:
                        mols[name] = mol
                except KeyboardInterrupt as err:
                    raise err
                except Exception as error:
                    if throw_on_error:
                        raise error
                    cls.journal.error(f'{error.__class__.__name__} for {name} - {error}')
        return mols

    @classmethod
    def extract_mol(cls,
                     name: str,
                     filepath: Optional[str] = None,
                     block: Optional[str] = None,
                     smiles: Optional[str] = None,
                     ligand_resn: str = 'LIG',
                     removeHs: bool = False,
                     proximityBonding: bool = False,
                     throw_on_error : bool = False) -> Chem.Mol:
        """
        Extracts the ligand of 3-name ``ligand_resn``
        from the PDB file ``filepath`` or from the PDB block ``block``.
        Corrects the bond order with SMILES if given.
        If there is a covalent bond with another residue the bond is kept as a ``*``/R.
        If the SMILES provided lacks the ``*`` element, the SMILES will be converted (if a warhead is matched),
        making the bond order correction okay.

        :param name: name of ligand
        :type name: str
        :param filepath: PDB file
        :type filepath: str
        :param smiles: SMILES
        :type smiles: str
        :param ligand_resn: 3letter PDB name of residue of ligand
        :type ligand_resn: str
        :param removeHs: Do you trust the hydrgens in the the PDB file?
        :type removeHs: bool
        :param throw_on_error: If an error occurs in the template step, raise error.
        :type throw_on_error: bool
        :return: rdkit Chem object
        :rtype: Chem.Mol
        """
        if filepath:
            readfun = Chem.MolFromPDBFile
            data = filepath
        elif block:
            readfun = Chem.MolFromPDBBlock
            data = block
        else:
            raise ValueError(f'Provide either a filepath or a block not neither')
        holo = readfun(data, proximityBonding=proximityBonding, removeHs=removeHs)
        if holo is None:
            cls.journal.warning(f'PDB {filepath} is problematic. Skipping sanitization.')
            holo = readfun(data, proximityBonding=False, removeHs=True, sanitize=False)
        mol = Chem.SplitMolByPDBResidues(holo, whiteList=[ligand_resn])[ligand_resn]
        mol = add_dummy_to_mol(mol, ligand_resn, holo)
        if smiles is not None:
            if '*' in Chem.MolToSmiles(mol) and '*' not in smiles:
                new_smiles = cls.make_covalent(smiles)
                if new_smiles:
                    cls.journal.info(f'{name} is covalent but ' +
                                      'a non covalent SMILES was passed, which was converted')
                    smiles = new_smiles
                else:
                    cls.journal.warning(f'{name} is covalent but ' +
                                        'a non covalent SMILES was passed, which failed to convert')
            else:
                pass
            try:
                template = Chem.MolFromSmiles(smiles)
                # template = AllChem.DeleteSubstructs(template, Chem.MolFromSmiles('*'))
                AllChem.SanitizeMol(mol)
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except ValueError as error:
                if throw_on_error:
                    raise error
                else:
                    cls.journal.warning(f'{name} failed at template-guided bond order correction ' +
                                        f'- ({type(error)}: {error}).')
        mol.SetProp('_Name', name)
        return mol


    # =================== From Files ===================================================================================

    @classmethod
    def from_files(cls, folder: str) -> _VictorUtils:  # future.annotation is active, so no cyclical.
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
        fragjson = os.path.join(folder, f'{self.long_name}.monster.json')
        if os.path.exists(fragjson):
            fd = json.load(open(fragjson))
            self.smiles = fd['smiles']
            self.is_covalent = True if '*' in self.smiles else False
            self.monster = Monster(self.hits)
            # self.monster.place(mol=self.mol,
            #                    attachment=None,
            #                    merging_mode='off')
            self.monster.positioned_mol = self.mol
            self.monster.positioned_mol.SetProp('_Origins', json.dumps(fd['origin']))
        else:
            self.is_covalent = None
            self.smiles = ''
            self.monster = None
            self.journal.info(f'{self.long_name} - no monster json')
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
        self.unminimized_pdbblock = None
        self.igor = None
        self.minimized_pdbblock = None
        # buffers etc.
        self._warned = []
        minjson = os.path.join(folder, f'{self.long_name}.minimised.json')
        self.mrmsd = mRMSD.mock()
        if os.path.exists(minjson):
            md = json.load(open(minjson))
            self.energy_score = md["Energy"]
            self.mrmsd.mrmsd = md["mRMSD"]
            self.mrmsd.rmsds = md["RMSDs"]
            self.igor = Igor.from_pdbfile(
                pdbfile=os.path.join(folder, self.long_name + '.holo_minimised.pdb'),
                params_file=os.path.join(folder, self.long_name + '.params'),
                constraint_file=os.path.join(folder, self.long_name + '.con'))
            # victor._fix_minimized adds to igor.mol_from_pose
            self.minimized_mol = self._fix_minimized()
        else:
            self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                                 'unbound_ref2015': {'total_score': float('nan')}}

            self.journal.info(f'{self.long_name} - no min json')
        minmol = os.path.join(folder, f'{self.long_name}.minimised.mol')
        if os.path.exists(minmol):
            self.minimized_mol = Chem.MolFromMolFile(minmol, sanitize=False, removeHs=False)
        else:
            self.minimized_mol = None
        return self

    # =================== Guess ===================================================================================

    @classmethod
    def get_isomers(cls, mol: Chem.Mol) -> List[Chem.Mol]:
        """
        For placement operations in particular it is important to differentiate the
        isomers. Therefore requiring multiple victor calls.
        """
        enumerate_tautomers = rdMolStandardize.TautomerEnumerator().Enumerate
        enumerate_stereoisomers = EnumerateStereoisomers.EnumerateStereoisomers
        return [tauto for stereo in enumerate_stereoisomers(mol)
                for tauto in enumerate_tautomers(stereo)
                ]

    @classmethod
    def get_isomers_smiles(cls, smiles: str) -> List[str]:
        """
        Same as `get_isomers`, but with smiles.
        """
        return list(map(Chem.MolToSmiles, cls.get_isomers(Chem.MolFromSmiles(smiles))))

    @classmethod
    def guess_warhead(cls, smiles: str) -> Tuple[str, str]:
        """
        Going backwards by guessing what the warhead is.
        Normally there'd be better data handling so no guessing
        """
        if not smiles or not isinstance(smiles, str):
            return '', 'noncovalent'
        if '*' not in smiles:
            return smiles, 'noncovalent'
        mol = Chem.MolFromSmiles(smiles)
        warhead_defs = []
        for warhead_def in cls.warhead_definitions:
            if mol.HasSubstructMatch(Chem.MolFromSmiles(warhead_def['covalent'])):
                warhead_defs.append(warhead_def)
        if not warhead_defs:
            cls.journal.warning(f'Could not match {smiles} to a warhead definition in `Victor.warhead_definitions`')
            return smiles.replace('*', 'O'), 'noncovalent'
        # most complex one first!
        warhead_def = sorted(warhead_defs, key=lambda d: -len(d['covalent_atomnames']))[0]
        unrxn_mols = AllChem.ReplaceSubstructs(mol=mol,
                                               query=Chem.MolFromSmiles(warhead_def['covalent']),
                                               replacement=Chem.MolFromSmiles(warhead_def['noncovalent'])
                                               )  # noqa it is filled.
        return Chem.MolToSmiles(unrxn_mols[0]), warhead_def['name']

    def get_plip_interactions(self):
        """
        Optional, but useful to have.
        And highly experimental!
        Get the interactions from PLIP.
        """
        from .plip import SerialPLIPper
        from plip.basic import config
        # PLIP is a bit problematic with hydrogens as it doesn't like them and loses atomtypes
        config.NOHYDRO = False  # default
        clean_block = ''
        for l in self.minimized_pdbblock.split('\n'):
            if l.startswith('#'):
                break
            if '    H' in l:
                continue
            clean_block += l + '\n'
        plipper = SerialPLIPper(pdb_block=clean_block,
                          resn=self.ligand_resn,
                          resi=int(self.ligand_resi[:-1]),
                          chain=self.ligand_resi[-1])
        setattr(self, 'plipper', plipper)  # new attribute
        interaction_set = plipper.get_interaction_set(plipper.pdb_block)
        return plipper.get_interaction_counts(interaction_set)

    @staticmethod
    def to_simple_smiles(mol: Chem.Mol) -> str:
        if not isinstance(mol, Chem.Mol) or not mol.GetNumAtoms():
            return ''
        try:
            mol = AllChem.RemoveAllHs(mol)
            for atom in mol.GetAtoms():
                atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            return Chem.MolToSmiles(mol)
        except:
            return ''

    def migrate_sw_origins(self, row: pd.Series) -> Dict[str, Dict[int, int]]:
        """
        Given a Victor object and a SmallWorld seach result row, return the "custom_map"
        """
        # hit names to [hit to minimised]
        origins: Dict[str, Dict[int, int]] = self.monster.convert_origins_to_custom_map(forbiddance=False)
        # minimised
        query2minimized: Dict[int, int] = dict(
            enumerate(self.minimized_mol.GetSubstructMatch(Chem.MolFromSmiles(row['qrySmiles']))))
        minimized2query: Dict[int, int] = dict(zip(query2minimized.values(), query2minimized.keys()))
        # candiate
        query2candidate: Dict[int, int] = dict(enumerate(row['atomMap']))
        candidate2query: Dict[int, int] = dict(zip(query2candidate.values(), query2candidate.keys()))
        # mixing
        minimized2candidate: Dict[int, int] = {m: query2candidate.get(q, -1) for m, q in minimized2query.items()}
        origins2candidate = {name: {o: minimized2candidate.get(m, -1) for o, m in origins[name].items()} for name in
                             origins}
        return origins2candidate

