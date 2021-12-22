from ._victor_igor import _VictorIgor

import os, re
from typing import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Constraints

# ----------------------------------------------------------------------------------------------------------------------

class _VictorCommon(_VictorIgor):

    def make_output_folder(self):
        path = os.path.join(self.work_path, self.long_name)
        if not os.path.exists(self.work_path):
            os.mkdir(self.work_path)
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            self.journal.warning(f'{self.long_name} - Folder {path} exists.')

    def _save_prerequisites(self):
        self._log_warnings()
        #  saving params
        self.journal.debug(f'{self.long_name} - saving params')
        params_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params')
        self.params.dump(params_file)
        # saving holo
        self.journal.debug(f'{self.long_name} - saving holo (unmimised)')
        holo_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_unminimised.pdb')
        with open(holo_file, 'w') as w:
            w.write(self.unminimized_pdbblock)
        # saving constraint
        if self.constraint is not None:
            self.journal.debug(f'{self.long_name} - saving constraint')
            constraint_file = os.path.join(self.work_path, self.long_name, self.long_name + '.con')
            self.constraint.dump(constraint_file)
        else:  # basically impossible.
            constraint_file = ''
        return params_file, holo_file, constraint_file

    # =================== Constraint & attachment ======================================================================

    def _get_constraint(self, extra_constraint: Optional[str] = None) -> Union[Constraints, None]:
        # deal with covalent and non covalent separately
        if self.is_covalent:
            self.journal.debug(f'{self.long_name} - is covalent.')
            constraint = self._fix_covalent()
            if extra_constraint:
                constraint.custom_constraint += self.extra_constraint
            return constraint
        else:
            self.journal.debug(f'{self.long_name} - is not covalent.')
            constraint = self._fix_uncovalent()
            if extra_constraint:
                constraint.custom_constraint += self.extra_constraint
            return constraint

    def _fix_uncovalent(self):
        return Constraints.mock()

    def _fix_covalent(self):
        self.journal.debug(f'{self.long_name} - fixing for covalent')
        # to make life easier for analysis, CX is the attachment atom, CY is the one before it.
        war_def = self._get_war_def()
        warhead = Chem.MolFromSmiles(war_def['covalent'])
        self.params.rename_by_substructure(warhead, war_def['covalent_atomnames'])
        cov_def = [d for d in self.covalent_definitions if d['residue'] == self.covalent_resn][0]
        self.journal.debug(f'{self.long_name} - has a {war_def["name"]}')
        cons = Constraints(smiles=(war_def['covalent'], cov_def['smiles']),
                           names=[*war_def['covalent_atomnames'], *cov_def['atomnames']],
                           ligand_res=self.ligand_resi,
                           target_res=self.covalent_resi)
        # user added constraint
        if 'constraint' in war_def:
            cons.custom_constraint = war_def['constraint']
        return cons

    def add_extra_constraint(self, new_constraint:Union[str]=None):
        if new_constraint is None:
            return # do nothing
        new_constraint = new_constraint.strip()
        if self.extra_constraint is None:
            self.extra_constraint = new_constraint
        self.extra_constraint = self.extra_constraint.strip() + '\n' + new_constraint.strip()

    def make_coordinate_constraints_for_placement(self,
                                    mol: Optional[Chem.Mol] = None,
                                    origins: Optional[List[List[str]]] = None,
                                    std: Optional[List[float]] = None,
                                    mx: Optional[List[float]] = None) -> str:
        """
        See also ``make_coordinate_constraints_for_combination`` in combine.
        This is the normal function and uses the origin data,
        while the other constrains based on lack of novel attribute.

        :param mol: self.monster.positioned_mol if ommitted
        :param origins: self.monster.origin_from_mol(self.monster.positioned_mol)  if ommitted
            list of list of names of hit atoms used as original position
        :param std: self.monster.stdev_from_mol(mol)(self.monster.positioned_mol)  if omitted
            list of standard devs
        :param mx: elf.monster.max_from_mol(self.monster.positioned_mol)  if omitted
            list of maximum euclidean distance
        :return:
        """
        lines = []
        if mol is None:
            mol = self.monster.positioned_mol
        if origins is None:
            origins = self.monster.origin_from_mol(mol)
        if std is None:
            std = self.monster.stdev_from_mol(mol)
        if mx is None:
            mx = self.monster.max_from_mol(mol)
        conf = self.monster.positioned_mol.GetConformer()
        # Calculate
        for i in range(mol.GetNumAtoms()):
            if len(origins[i]) > 0:
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == '*':
                    continue
                elif atom.GetPDBResidueInfo() is None:
                    self.journal.critical(f'Atom {i} ({atom.GetSymbol()}) has no name!')
                    continue
                pos = conf.GetAtomPosition(i)
                if self.constraint_function_type.upper() == 'HARMONIC':
                    fxn = f'HARMONIC 0 {std[i] + 1}'
                elif self.constraint_function_type.upper() == 'FLAT_HARMONIC':
                    if len(origins[i]) > 1:
                        fxn = f'FLAT_HARMONIC 0 1.0 {mx[i]}'
                    else:
                        fxn = f'HARMONIC 0 1.0'
                elif self.constraint_function_type.upper() == 'BOUNDED':
                    fxn = f'BOUNDED 0 {mx[i]} 1 0.5 TAG'
                else:
                    raise ValueError(f'{self.constraint_function_type} is not HARMONIC or FADE or BOUNDED')
                atomname = atom.GetPDBResidueInfo().GetName()
                lines.append(f'CoordinateConstraint {atomname} {self.ligand_resi} ' +
                             f'CA {self.covalent_resi} ' +
                             f'{pos.x} {pos.y} {pos.z} {fxn}\n')
        return ''.join(lines)

    def make_coordinate_constraints_for_combination(self):
        """
        See also ``cls.make_coordinate_constraints_for_placement``.
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

    # ------------------------------------------------------------------------------------------------------------------

    def _get_attachment_from_pdbblock(self) -> Union[None, Chem.Mol]:
        """
        Yes, yes, I see the madness in using pymol to get an atom for rdkit to make a pose for pyrosetta.
        Hence why `find_attachment` will replace it.
        todo `_get_attachment_from_pdbblock` --> `find_attachment`
        """
        import pymol2

        self.journal.debug(f'{self.long_name} - getting attachemnt atom')
        if not self.covalent_resn:
            return None
        else:
            if isinstance(self.covalent_resi, str):
                resi, chain = re.match('(\d+)(\w)', self.covalent_resi).groups()
                resi = int(resi)
            else:
                resi = self.covalent_resi
                chain = None
            with pymol2.PyMOL() as pymol:
                pymol.cmd.read_pdbstr(self.apo_pdbblock, 'prot')
                if self.covalent_resn == 'CYS':
                    name = 'SG'
                else:
                    raise NotImplementedError('only done for cys atm')
                try:
                    if chain is not None:
                        pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name} and chain {chain}')
                    else:
                        pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name}')
                except:
                    pdb = pymol.cmd.get_pdbstr(f'resi {resi} and name {name}')
                return Chem.MolFromPDBBlock(pdb)

    def _get_war_def(self):
        for war_def in self.warhead_definitions:
            warhead = Chem.MolFromSmiles(war_def['covalent'])
            if self.mol.HasSubstructMatch(warhead):
                return war_def
        else:
            if self.mol.HasSubstructMatch(Chem.MolFromSmiles('*C')):
                self.journal.warning('Unknown type of covalent')
                return {'name': 'unknown',
                        'covalent': 'C~C*',
                        'covalent_atomnames': ['CY', 'CX', 'CONN1'],
                        'noncovalent': 'C~[C+]',  # clearly not
                        'noncovalent_atomnames': ['CY', 'CX']
                        }
            else:
                raise ValueError(f'{self.long_name} - Unsure what the warhead is {self.smiles}.')

    @classmethod
    def inventorize_warheads(cls, hits: List[Chem.Mol], covalent_form: bool = True) -> List[str]:
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
    def _get_warhead_from_name(cls, warhead_name: str, covalent_form: bool) -> Chem.Mol:
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
    def _get_warhead_from_definition(cls, war_def: dict, covalent_form: bool):
        if covalent_form:
            wh = Chem.MolFromSmiles(war_def['covalent'])
        else:
            wh = Chem.MolFromSmiles(war_def['noncovalent'])
        return wh

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

    def harmonize_warheads(self, hits, warhead_harmonisation, covalent_form=True):
        """
        Harmonises and marks the atoms with `_Warhead` Prop.

        :param hits:
        :param warhead_harmonisation:
        :param covalent_form:
        :return:
        """
        inventory = self.inventorize_warheads(hits, covalent_form)
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
