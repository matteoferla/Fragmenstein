from ._victor_base_mixin import _VictorBaseMixin
from ..monster import Monster
from ..igor import Igor
from ..m_rmsd import mRSMD
from ..rectifier import Rectifier
from typing import List, Optional, Dict, Union, Callable
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
import time, warnings


class _VictorAutomergeMixin(_VictorBaseMixin):


    @classmethod
    def combine(cls,
                hits: List[Chem.Mol],
                pdb_filename: str,
                ligand_resn: str = 'LIG',
                ligand_resi: Union[int, str] = '1B',
                covalent_resn: str = 'CYS',  # no other option is accepted.
                covalent_resi: Optional[Union[int, str]] = None,
                extra_constraint: Union[str] = None,
                pose_fx: Optional[Callable] = None,
                atomnames: Optional[Dict[int, str]] = None,
                warhead_harmonisation: str='first'
                ):
        """
        Combines the hits without a template.
        If the class attribute ``monster_throw_on_discard`` is True, it will raise an exception if it cannot.

        The cutoff distance is controlled by class attribute ``monster_joining_cutoff``.
        At present this just adds a hydrocarbon chain, no fancy checking for planarity.

        The hits are collapsed, merged, expanded and bonded by proximity.
        In ``(self.monster.expand_ring(..., bonded_as_original=False)`` changing to True, might work, but most likely won't.

        ``warhead_harmonisation`` fixes the warhead in the hits to be homogeneous.

        * ``keep``. Don't do anything
        * ``none``. strip warheads
        * ``first``. Use first warhead
        * warhead name. Use this warhead.

        :param hits:
        :param pdb_filename:
        :param ligand_resn:
        :param ligand_resi:
        :param covalent_resn:
        :param covalent_resi:
        :param extra_constraint:
        :param pose_fx:
        :param atomnames:
        :param warhead_harmonisation: keep | strip | first | chloracetimide | nitrile ...
        :return:
        """
        self = cls.__new__(cls)
        self.monster_merging_mode = 'full' # needed solely for logkeeping
        self.long_name = '-'.join([h.GetProp('_Name') for h in hits])
        with open(pdb_filename) as fh:
            self.apo_pdbblock = fh.read()
        self.journal.debug(f'{self.long_name} - harmonising warheads on hits in "{warhead_harmonisation}" mode')
        with warnings.catch_warnings(record=True) as self._warned:
            self.hits = self.harmonise_warheads(hits, warhead_harmonisation, covalent_form=True)
            self._log_warnings()
        self.ligand_resn = ligand_resn.upper()
        self.ligand_resi = ligand_resi
        self.covalent_resn = covalent_resn.upper()
        self.covalent_resi = covalent_resi
        self.atomnames = atomnames
        self.extra_constraint = extra_constraint
        self.pose_fx = pose_fx
        # these are calculated
        starhits = any(['*' in Chem.MolToSmiles(h) for h in hits])
        if starhits and (self.covalent_resi is None or self.covalent_resn is None):
            raise ValueError(f'{self.long_name} - is covalent but without known covalent residues')
        elif warhead_harmonisation == 'strip':
            self.is_covalent = False
        elif starhits:
            self.is_covalent = True
        else:
            self.is_covalent = False
        self.params = None
        self.mol = None
        self.smiles = None
        self.constraint = None
        self.monster = None
        # ====== debug: absent in main mode.
        self.modifications = []  # list of the various steps during fragment merging mode.
        # ======
        self.unminimised_pdbblock = None
        self.igor = None
        self.minimised_pdbblock = None
        self.minimised_mol = None
        if self.monster_average_position:
            # I need to code this case.
            self.journal.warning('`monster_average_position` class attribute == True does not apply here')
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.tick = time.time()
        self.tock = float('inf')
        self.error = ''
        self._safely_do(execute=self._combine_main,
                        resolve=self._resolve,
                        reject=self._reject)
        return self

    def _combine_main(self):
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self.monster = Monster(hits=self.hits,
                               debug_draw=self.monster_debug_draw,
                               average_position=self.monster_average_position
                               )
        self.monster.modifications = self.modifications
        self.monster.merge(keep_all=self.monster_throw_on_discard,
                           collapse_rings=True,
                           joining_cutoff=self.monster_joining_cutoff # Å
                            )
        self.mol = self.monster.positioned_mol
        self.smiles = Chem.MolToSmiles(self.mol)
        if self.monster_debug_draw:
            picture = Chem.CombineMols(Chem.CombineMols(self.hits[0], self.hits[1]), self.monster.positioned_mol)
            AllChem.Compute2DCoords(picture)
            self.monster.draw_nicely(picture)
        # making folder.
        self._make_output_folder()
        # paramterise
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        self.params = Params.load_mol(self.mol, name=self.ligand_resn)
        self.params.NAME = self.ligand_resn # force it.
        self.params.polish_mol()
        # get constraint
        self.constraint = self._get_constraint(self.extra_constraint)
        self.constraint.custom_constraint += self._make_coordinate_constraints_for_unnovels()
        # _get_constraint will have changed the names in params.mol so the others need changing too!
        # namely  self.params.rename_by_substructure happend.
        self.mol = Chem.Mol(self.params.mol)
        self.monster.positioned_mol = Chem.Mol(self.mol)
        # those Hs lack correct names and charge!!
        self.params.add_Hs()
        self.params.convert_mol()
        # self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        # self.params.CHI.data = []  # TODO check if chi fix is okay
        self._log_warnings()
        self.post_params_step()
        self.monster_merging_mode = 'full'
        self.unminimised_pdbblock = self._plonk_monster_in_structure()
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.unbound_pose = self.params.test()
        self._checkpoint_alpha()
        self._checkpoint_bravo()
        self.igor = Igor.from_pdbblock(pdbblock=self.unminimised_pdbblock,
                                       params_file=params_file,
                                       constraint_file=constraint_file,
                                       ligand_residue=self.ligand_resi,
                                       key_residues=[self.covalent_resi])
        # user custom code.
        if self.pose_fx is not None:
            self.journal.debug(f'{self.long_name} - running custom pose mod.')
            self.pose_fx(self.igor.pose)
        else:
            self.pose_mod_step()
        # storing a roundtrip
        self.unminimised_pdbblock = self.igor.pose2str()
        # minimise until the ddG is negative.
        self.reanimate_n_store()
        self.journal.debug(f'{self.long_name} - Completed')

    def reanimate_n_store(self):
        ddG = self.reanimate()
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {ddG} kcal/mol {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        self.journal.debug(f'{self.long_name} - Completed')

    def _make_coordinate_constraints_for_unnovels(self):
        """
        See also ``cls._make_coordinate_constraints``.
        This operates based on ``atom.HasProp('_Novel')``, not origins!
        :return:
        """
        lines = []
        conf = self.monster.positioned_mol.GetConformer()
        for i, atom in enumerate(self.monster.positioned_mol.GetAtoms()):
            if atom.GetSymbol() == '*':
                continue
            elif atom.HasProp('_Novel') and atom.GetBoolProp('_Novel'):
                continue # novels
            elif atom.GetPDBResidueInfo() is None:
                self.journal.critical(f'Atom {i} ({atom.GetSymbol()}) has no name!')
                continue
            else:
                pos = conf.GetAtomPosition(i)
                fxn = f'HARMONIC 0 1' # the other do not make sense here.
                lines.append(f'CoordinateConstraint {atom.GetPDBResidueInfo().GetName()} {self.ligand_resi} ' + \
                             f'CA {self.covalent_resi} ' + \
                             f'{pos.x} {pos.y} {pos.z} {fxn}\n')
        return ''.join(lines)

    @classmethod
    def inventorise_warheads(cls, hits: List[Chem.Mol], covalent_form: bool=True) -> List[str]:
        """
        Get the warhead types of the list of hits

        :param hits:
        :param covalent_form: Are the hits already covalent (with *)
        :return: list of non-covalent, chloroacetimide, etc.
        """
        inventory = ['noncovalent'] * len(hits)
        for war_def in cls.warhead_definitions:
            wh = cls._get_warhead_from_definition(war_def, covalent_form)
            for i, hit in enumerate(hits):
                if hit.HasSubstructMatch(wh):
                    inventory[i] = war_def['name']
        return inventory

    @classmethod
    def _get_warhead_from_name(cls, warhead_name:str, covalent_form:bool) -> Chem.Mol:
        """
        get_warhead_definition returns a definition, this retursn a mol.

        :param warhead_name:
        :param covalent_form:
        :return:
        """
        war_def = cls.get_warhead_definition(warhead_name)
        wh = cls._get_warhead_from_definition(war_def, covalent_form)
        return wh

    @classmethod
    def _get_warhead_from_definition(cls, war_def:dict, covalent_form:bool):
        if covalent_form:
            wh = Chem.MolFromSmiles(war_def['covalent'])
        else:
            wh = Chem.MolFromSmiles(war_def['noncovalent'])
        return wh

    def harmonise_warheads(self, hits, warhead_harmonisation, covalent_form=True):
        """
        Harmonises and marks the atoms with `_Warhead` Prop.

        :param hits:
        :param warhead_harmonisation:
        :param covalent_form:
        :return:
        """
        inventory = self.inventorise_warheads(hits, covalent_form)
        # mark warhead atoms.
        for hit, warhead_name in zip(hits, inventory):
            if warhead_name != 'noncovalent':
                wh = self._get_warhead_from_name(warhead_name, covalent_form)
                match = hit.GetSubstructMatch(wh)
                if match == ():
                    self.journal.warning(f'{self.long_name} - failed to mark warhead. What is it??')
                else:
                    for i in match:
                        hit.GetAtomWithIdx(i).SetBoolProp('_Warhead', True)
        # harmonise
        if warhead_harmonisation == 'keep':
            return hits
        elif warhead_harmonisation == 'strip':
            new_hits = []
            for hit, warhead_name in zip(hits, inventory):
                if warhead_name != 'noncovalent':
                    wh = self._get_warhead_from_name(warhead_name, covalent_form)
                    nhit = AllChem.DeleteSubstructs(hit, wh)
                    Chem.SanitizeMol(nhit)
                    new_hits.append(nhit)
                else:
                    new_hits.append(hit)

            return new_hits
        elif warhead_harmonisation == 'first':
            if len(set(inventory) - {'noncovalent'}) <= 1:
                return hits
            else:
                first = None
                new_hits = []
                for hit, warhead_name in zip(hits, inventory):
                    if warhead_name != 'noncovalent':
                        if first is None:
                            first = warhead_name
                            new_hits.append(hit)
                        elif warhead_name == first:
                            new_hits.append(hit)
                        else:
                            wh = self._get_warhead_from_name(warhead_name, covalent_form)
                            nhit = AllChem.DeleteSubstructs(hit, wh)
                            Chem.SanitizeMol(nhit)
                            new_hits.append(nhit)
                    else:
                        new_hits.append(hit)
                return new_hits
        else: # it is a warhead name.
            raise NotImplementedError




