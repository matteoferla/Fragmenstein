########################################################################################################################
__doc__ = \
    """
These are extras for the Monster step
    """

########################################################################################################################

from typing import List, Optional, Tuple, Dict
from warnings import warn

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, Draw

import json

try:
    from IPython.display import SVG, display
except (KeyError, ImportError):
    warn('No Jupyter notebook installed. `.draw_nicely` will not work.')
    SVG = lambda *args, **kwargs: print('Install IPython...')
    display = lambda *args, **kwargs: print('Install IPython...')

from ._communal import _MonsterCommunal
from .positional_mapping import GPM


########################################################################################################################


class _MonsterUtil(_MonsterCommunal, GPM):

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
            followup_placed = cls.positioned_mol
        if hits is None:  # instance
            assert hasattr(cls, '__class__'), 'if called as a classmethod the list of hits need to be provided.'
            hits = cls.hits
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
            elif atom.HasProp('_ori_name'): # single name.
                origin.append([atom.GetProp('_ori_name')])
            else:
                origin.append([])
        return origin

    def guess_origins(self, mol: Chem.Mol = None, hits: Optional[List[Chem.Mol]] = None):
        """
        Given a positioned mol guess its origins...

        :param mol:
        :return:
        """

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

    def make_pse(self, filename='test.pse', extra_mols: Optional[Chem.Mol] = None):
        """
        This is specifically for debugging the full fragment merging mode.
        For general use. Please use the Victor method ``make_pse``.

        :param filename:
        :return:
        """
        assert '.pse' in filename, 'Must be a pymol pse extension!'
        import pymol2
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
            if extra_mols:
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
            print(mol.GetProp('_Name'))
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
        except Exception as err:
            warn(f'*{err.__class__.__name__}* : {err}')
            display(x)

    def mmff_minimize(self, mol: Optional[Chem.Mol] = None) -> None:
        """
        Minimises a mol, or self.positioned_mol if not provided, with MMFF constrained to 2 Å.
        Gets called by Victor if the flag .monster_mmff_minimisation is true during PDB template construction.

        :param mol: opt. mol. modified in place.
        :return: None
        """
        success = True
        if mol is None and self.positioned_mol is None:
            raise ValueError('No valid molecule')
        elif mol is None:
            mol = self.positioned_mol
        else:
            pass  # mol is fine
        # protect
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0)):
            atom.SetBoolProp('_IsDummy', True)
            atom.SetAtomicNum(16)
        #
        mol.UpdatePropertyCache()
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
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_Novel', negate=True)):
            i = atom.GetIdx()
            ff.MMFFAddPositionConstraint(i, 2, 10)
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_IsDummy')):
            i = atom.GetIdx()
            ff.MMFFAddPositionConstraint(i, 0.1, 10)
        try:
            m = ff.Minimize()
            if m == -1:
                self.journal.error('MMFF Minisation could not be started')
            elif m == 0:
                self.journal.info('MMFF Minisation was successful')
            elif m == 1:
                self.journal.info('MMFF Minisation was run, but the minimisation was not unsuccessful')
            else:
                self.journal.critical("Iä! Iä! Cthulhu fhtagn! Ph'nglui mglw'nafh Cthulhu R'lyeh wgah'nagl fhtagn")
        except RuntimeError as error:
            self.journal.error(f'MMFF minimisation failed {error.__class__.__name__}: {error}')
            success = False
        # deprotect
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_IsDummy')):
            atom.SetAtomicNum(0)
        return success
