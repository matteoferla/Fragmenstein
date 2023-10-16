########################################################################################################################
__doc__ = \
    """
These are extras for the Monster step
    """

########################################################################################################################
from itertools import combinations
from typing import List, Optional, Tuple, Dict, Union, Iterable
from warnings import warn

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, Draw, rdMolAlign, rdqueries

import json

try:
    from IPython.display import SVG, display
except (KeyError, ImportError):
    warn('No Jupyter notebook installed. `.draw_nicely` will not work.')
    SVG = lambda *args, **kwargs: print('Install IPython...')
    display = lambda *args, **kwargs: print('Install IPython...')

from ._communal import _MonsterCommunal
from .positional_mapping import GPM
from ._util_compare import _MonsterUtilCompare


########################################################################################################################


class _MonsterUtil(_MonsterCommunal, GPM, _MonsterUtilCompare):

    @classmethod
    def get_combined_rmsd(cls, followup_moved: Chem.Mol, followup_placed: Optional[Chem.Mol] = None,
                          hits: Optional[List[Chem.Mol]] = None) -> float:
        """
        Depracated.
        The inbuilt RMSD calculations in RDKit align the two molecules, this does not align them.
        This deals with the case of multiple hits.
        For euclidean distance the square root of the sum of the differences in each coordinates is taken.
        For a regular RMSD the still-squared distance is averaged before taking the root.
        Here the average is done across all the atom pairs between each hit and the followup.
        Therefore, atoms in followup that derive in the blended molecule by multiple atom are scored multiple times.

        As a classmethod ``followup_placed`` and ``hits`` must be provided. But as an instance method they don't.

        :param followup_moved: followup compound moved by Igor or similar
        :param followup_placed: followup compound as placed by Monster
        :param hits: list of hits.
        :return: combined RMSD
        """
        # class or instance?
        if followup_placed is None:  # instance
            assert hasattr(cls, '__class__'), 'if called as a classmethod the list of hits need to be provided.'
            followup_placed = cls.positioned_mol  # noqa it's not a class but an instance
        if hits is None:  # instance
            assert hasattr(cls, '__class__'), 'if called as a classmethod the list of hits need to be provided.'
            hits = cls.hits  # noqa it's not a class but an instance
        for i in range(followup_placed.GetNumAtoms()):
            assert followup_placed.GetAtomWithIdx(i).GetSymbol() == followup_moved.GetAtomWithIdx(
                i).GetSymbol(), 'The atoms order is changed.'
        if followup_moved.GetNumAtoms() > followup_placed.GetNumAtoms():
            warn(
                f'Followup moved {followup_moved.GetNumAtoms()} has more atoms that followup placed {followup_placed.GetNumAtoms()}. Assuming these are hydrogens.')
        # calculate
        tatoms = 0
        d = 0
        for hit in hits:
            mapping = list(cls.get_positional_mapping(followup_placed, hit).items())
            tatoms += len(mapping)
            if len(mapping) == 0:
                continue
            d += cls._get_square_deviation(followup_moved, hit, mapping)
        return d / tatoms ** 0.5

    @classmethod
    def get_pair_rmsd(cls, molA, molB, mapping: List[Tuple[int, int]]) -> float:
        return (cls._get_square_deviation(molA, molB, mapping) / len(mapping)) ** 0.5

    def _get_square_deviation(self, molA, molB, mapping):
        confA = molA.GetConformer()
        confB = molB.GetConformer()
        return sum([(confA.GetAtomPosition(a).x - confB.GetAtomPosition(b).x) ** 2 +
                    (confA.GetAtomPosition(a).y - confB.GetAtomPosition(b).y) ** 2 +
                    (confA.GetAtomPosition(a).z - confB.GetAtomPosition(b).z) ** 2 for a, b in mapping])

    @property
    def num_common(self) -> int:
        template = self._get_last_template()
        mcs = rdFMCS.FindMCS([template, self.initial_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder)
        return Chem.MolFromSmarts(mcs.smartsString).GetNumAtoms()

    def _get_last_template(self):
        if 'chimera' in self.modifications:
            template = self.modifications['chimera']
        elif 'scaffold' in self.modifications:
            template = self.modifications['scaffold']
        else:
            raise KeyError('There is no chimeric or reg scaffold/template to compare to.')

    @property
    def percent_common(self) -> int:
        return round(self.num_common / self.initial_mol.GetNumAtoms() * 100)

    def stdev_from_mol(self, mol: Chem.Mol = None):
        """
        these values are stored from Monster for scaffold, chimera and positioned_mol

        :param mol: Chem.Mol
        :return: stdev list for each atom
        """
        if mol is None:
            mol = self.positioned_mol
        return [atom.GetDoubleProp('_Stdev') if atom.HasProp('_Stdev') else 0 for atom in mol.GetAtoms()]

    def max_from_mol(self, mol: Chem.Mol = None):
        if mol is None:
            mol = self.positioned_mol
        return [atom.GetDoubleProp('_Max') if atom.HasProp('_Max') else 0 for atom in mol.GetAtoms()]

    def origin_from_mol(self, mol: Chem.Mol = None):
        """
        these values are stored from Monster for scaffold, chimera and positioned_mol
        See `make_chimera` or `place_from_map` for more info on _Origin

        :param mol: Chem.Mol
        :return: stdev list for each atom
        """
        if mol is None:
            mol = self.positioned_mol
        if mol.HasProp('_Origins'):
            return json.loads(mol.GetProp('_Origins'))
        origin = []
        for atom in mol.GetAtoms():
            # {'__computedProps': <rdkit.rdBase....>,
            # '_ori_i': 100, '_ori_name': 'x10976', '_x': 13.792, '_y': -1.785, '_z': 21.231,
            # '_GasteigerCharge': 0.0, '_GasteigerHCharge': 0.0, '_rType': 'Cl'}
            if atom.HasProp('_Origin'):
                x = atom.GetProp('_Origin')
                if x == 'none':
                    origin.append([])
                else:
                    origin.append(json.loads(x))
            elif atom.HasProp('_ori_name'):  # single name.
                origin.append([atom.GetProp('_ori_name') + '.' + atom.GetProp('_ori_i')])
            else:
                origin.append([])
        return origin

    def guess_origins(self, mol: Chem.Mol = None, hits: Optional[List[Chem.Mol]] = None):
        """
        Given a positioned mol guess its origins...

        :param mol:
        :return:
        """
        self.journal.debug('Guess origins called')
        if hits is None:
            hits = self.hits
        mappings = []
        for h, hit in enumerate(hits):
            hname = hit.GetProp('_Name')
            for hi, mi in self.get_positional_mapping(hit, mol).items():
                atom = mol.GetAtomWithIdx(mi)
                if atom.HasProp('_Novel') and atom.GetBoolProp('_Novel') == True:
                    continue  # flagged to avoid.
                elif atom.HasProp('_Origin') and atom.GetProp('_Origin') != 'none':
                    origin = json.loads(atom.GetProp('_Origin'))
                else:
                    origin = []
                origin.append(f'{hname}.{hi}')
                atom.SetProp('_Origin', json.dumps(origin))

    # class attribute for next method
    _i = 0

    def save_temp(self, mol):
        """
        This is a silly debug-by-print debug method. drop it in where you want to spy on stuff.
        """
        Chem.MolToMolFile(mol, f'debug_temp{self.i}.mol', kekulize=False)
        self._i += 1

    def save_commonality(self, filename: Optional[str] = None):
        """
        Saves an SVG of the followup fragmenstein monster with the common atoms with the chimeric scaffold highlighted.

        :param filename: optinal filename to save it as. Otherwise returns a Draw.MolDraw2DSVG object.
        :return:
        """
        template = self._get_last_template()
        mcs = rdFMCS.FindMCS([template, self.positioned_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder,
                             ringMatchesRingOnly=True)
        common = Chem.MolFromSmarts(mcs.smartsString)
        match = self.positioned_mol.GetSubstructMatch(common)
        d = self.draw_nicely(self.positioned_mol, show=False, highlightAtoms=match)
        if filename is None:
            return d
        else:
            assert '.svg' in filename, 'Can only save SVGs.'
            with open(filename, 'w') as w:
                w.write(d.GetDrawingText())

    def make_pse(self, filename='test.pse', extra_mols: Optional[Iterable[Chem.Mol]] = None):
        """
        This is specifically for debugging the full fragment merging mode.
        For general use. Please use the Victor method ``make_pse``.

        :param filename:
        :return:
        """
        assert '.pse' in filename, 'Must be a pymol pse extension!'
        import pymol2  # noqa it's an optional
        with pymol2.PyMOL() as pymol:
            tints = iter(
                ['wheat', 'palegreen', 'lightblue', 'paleyellow', 'lightpink', 'palecyan', 'lightorange', 'bluewhite'])
            # pymol.cmd.bg_color('white')
            for h, hit in enumerate(self.hits):
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit, kekulize=False), f'hit{h}')
                pymol.cmd.color(next(tints), f'hit{h} and name C*')
            if 'scaffold' in self.modifications:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.modifications['scaffold'], kekulize=False), f'scaffold')
                pymol.cmd.color('tv_blue', f'scaffold and name C*')
            if 'chimera' in self.modifications:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.modifications['chimera'], kekulize=False), f'chimera')
                pymol.cmd.color('cyan', f'chimera and name C*')
            if self.positioned_mol:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.positioned_mol, kekulize=False), f'followup')
                pymol.cmd.color('tv_green', f'followup and name C*')
            if self.mol_options:
                for i, mol in enumerate(self.mol_options):
                    pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), f'opt{i}')
                    pymol.cmd.color('grey50', f'opt{i} and name C*')
            pymol.cmd.hide('sticks')
            pymol.cmd.hide('cartoon')  # there should not be....
            pymol.cmd.show('lines', 'not polymer')
            if 'chimera' in self.modifications:
                pymol.cmd.show('sticks', 'chimera')
            if self.positioned_mol:
                pymol.cmd.show('sticks', 'followup')
            if extra_mols is not None:
                for mol in extra_mols:
                    name = mol.GetProp('_Name')
                    pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), name)
                    pymol.cmd.color('magenta', f'{name} and name C*')
            pymol.cmd.save(filename)

    def draw_nicely(self, mol, show=True, **kwargs) -> Draw.MolDraw2DSVG:
        """
        Draw with atom indices for Jupyter notebooks.


        :param mol:
        :param kwargs: Key value pairs get fed into ``PrepareAndDrawMolecule``.
        :return:
        """
        if mol.HasProp('_Name'):
            print(mol.GetProp('_Name'))  # this is a legit print not a rogue one
        d = Draw.MolDraw2DSVG(400, 400)
        d.drawOptions().addAtomIndices = True
        d.drawOptions().addStereoAnnotation = True
        d.drawOptions().prepareMolsBeforeDrawing = False
        d.drawOptions().dummiesAreAttachments = True
        x = Chem.Mol(mol)
        AllChem.Compute2DCoords(x)
        Chem.SanitizeMol(x, catchErrors=True)
        try:
            # x = Chem.MolFromSmiles(Chem.MolToSmiles(x, kekuleSmiles=False), sanitize=False)
            Draw.PrepareAndDrawMolecule(d, x, **kwargs)
            d.FinishDrawing()
            if show:
                display(SVG(d.GetDrawingText()))
            return d
        except KeyboardInterrupt as err:
            raise err
        except Exception as err:
            warn(f'*{err.__class__.__name__}* : {err}')
            display(x)

    def mmff_minimize(self, mol: Optional[Chem.Mol] = None,
                      ff_dist_thr: float = 5.,
                      ff_constraint: int = 10,
                      allow_lax: bool = True) -> bool:
        """
        Minimises a mol, or self.positioned_mol if not provided, with MMFF constrained to ff_dist_thr Å.
        Gets called by Victor if the flag .monster_mmff_minimisation is true during PDB template construction.

        :param mol: opt. mol. modified in place.
        :param ff_dist_thr: Distance threshold (Å) for atomic positions mapped to hits for  MMFF constrains.
                            if NaN then fixed point constraints (no movement) are used.
        :param ff_constraint: Force constant for MMFF constraints.
        :allow_lax: If True and the minimisation fails, the constraints are halved and the minimisation is rerun.
        :return: None

        Note that most methods calling this via Victor
        now use its ``.settings['ff_dist_thr']`` and ``.settings['ff_constraint']``
        and do not use the defaults.
        """
        success = True
        if mol is None and self.positioned_mol is None:
            raise ValueError('No valid molecule')
        elif mol is None:
            mol = self.positioned_mol
        else:
            pass  # mol is fine
        fixed_mode = str(ff_dist_thr).lower() == 'nan'
        # store for later (drift prevention)
        original_mol = Chem.Mol(mol)
        # protect (DummyMasker could be used here)
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0)):
            atom.SetBoolProp('_IsDummy', True)
            atom.SetAtomicNum(16)
        #
        # Chem.GetSymmSSSR(mol)
        # Chem.MolToMolFile(mol, 'test.mol')
        Chem.SanitizeMol(mol)
        #
        p = AllChem.MMFFGetMoleculeProperties(mol, 'MMFF94')
        if p is None:
            self.journal.error(f'MMFF cannot work on a molecule that has errors!')
            return False
        ff = AllChem.MMFFGetMoleculeForceField(mol, p)
        # restrain
        restrained = []
        novels = list(mol.GetAtomsMatchingQuery(rdqueries.HasPropQueryAtom('_Novel', negate=True)))
        if len(novels) == 0 and fixed_mode:
            self.journal.warning('No novel atoms found in fixed_mode (ff_dist_thr == NaN), ' + \
                                 'this is probably a mistake')
            return True  # nothing to do
        for atom in novels:
            i = atom.GetIdx()
            if atom.GetAtomicNum() == 1:
                continue
            if fixed_mode:
                ff.AddFixedPoint(i)
            else:
                ff.MMFFAddPositionConstraint(i, ff_dist_thr, ff_constraint)
            restrained.append(i)
        for atom in mol.GetAtomsMatchingQuery(rdqueries.HasPropQueryAtom('_IsDummy')):
            i = atom.GetIdx()
            ff.MMFFAddPositionConstraint(i, 0.1, ff_constraint)
            restrained.append(i)
        try:
            m = ff.Minimize()
            if m == -1:
                self.journal.error('MMFF Minisation could not be started')
                success = False
            elif m == 0:
                self.journal.info('MMFF Minisation was successful')
                success = True
            elif m == 1:
                self.journal.info('MMFF Minisation was run, but the minimisation was not unsuccessful')
                success = False
            else:
                self.journal.critical("Iä! Iä! Cthulhu fhtagn! Ph'nglui mglw'nafh Cthulhu R'lyeh wgah'nagl fhtagn")
                success = False
        except RuntimeError as error:
            self.journal.info(f'MMFF minimisation failed {error.__class__.__name__}: {error}')
            success = False
        if fixed_mode and len(restrained) > 0:
            # no need to align nothing could have moved
            return success
        if not success and allow_lax:
            success: bool = self.mmff_minimize(mol,
                                               ff_dist_thr=ff_dist_thr / 2,
                                               ff_constraint=ff_constraint // 2,
                                               allow_lax=False)
        # deprotect
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_IsDummy')):
            atom.SetAtomicNum(0)
        # prevent drift:
        rdMolAlign.AlignMol(mol, original_mol, atomMap=list(zip(restrained, restrained)))
        return success

    def MMFF_score(self, mol: Optional[Chem.Mol] = None, delta: bool = False) -> float:
        """
        Merck force field. Chosen over Universal for no reason at all.

        :param mol: ligand
        :type mol: Chem.Mol optional. If absent extracts from pose.
        :param delta: report difference from unbound (minimized)
        :type delta: bool
        :return: kcal/mol
        :rtype: float

        :warning: This was moved out of Igor. Victor has the method for calling it with igor.mol_from_pose
        """
        if mol is None:
            mol = self.positioned_mol
        try:
            mol = Chem.Mol(mol)  # copy!
            AllChem.UFFGetMoleculeForceField(mol)
            ff = AllChem.UFFGetMoleculeForceField(mol)
            ff.Initialize()
            # print(f'MMFF: {ff.CalcEnergy()} kcal/mol')
            if delta:
                pre = ff.CalcEnergy()
                ff.Minimize()
                post = ff.CalcEnergy()
                return pre - post
            else:
                return ff.CalcEnergy()
        except RuntimeError as err:
            warn(f'{err.__class__.__name__}: {err} (It is generally due to bad sanitisation)')
            return float('nan')

    def _get_substructure_from_idxs(self, mol: Chem.Mol, atomIdx_list: List[int]) -> \
            Tuple[Union[Chem.Mol, Chem.Mol], Dict[int, int]]:
        '''
        Given a molecule, extract the substructure molecule given selected atom idxs.
        :param mol:
        :param atomIdx_list:
        :return:
        '''
        bonds = []
        atommap = {}

        for i, j in combinations(atomIdx_list, 2):
            b = mol.GetBondBetweenAtoms(i, j)
            if b:
                bonds.append(b.GetIdx())

        newMol = Chem.PathToSubmol(mol, bonds, atomMap=atommap)
        if len(atommap) == 0:
            newMol = None
        return newMol, atommap

    @staticmethod
    def renumber_followup_custom_map(original_mol: Chem.Mol,
                                     new_mol: Chem.Mol,
                                     custom_map: Dict[str, Dict[int, int]]) -> Dict[str, Dict[int, int]]:
        """
        Give a followup ``original_mol`` and a copy but with its atom indices changed
        return a new map for the copy
        """
        # GetSubstructMatch returns the indices of the caller molecule, not the query molecule
        origin2new = {o: n for n, o in enumerate(original_mol.GetSubstructMatch(new_mol))}
        new_custom_map: Dict[str, Dict[int, int]] = {}
        for hit_name in custom_map:
            # replace the values of the dict with the new indices or with itself if negative (a forbidden atom)
            new_custom_map[hit_name] = {k: origin2new[v] if v >= 0 else v for k, v in custom_map[hit_name].items()}
        return new_custom_map
