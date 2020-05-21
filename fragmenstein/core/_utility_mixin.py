from typing import List, Optional, Tuple
from warnings import warn

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, Draw

import json

try:
    from IPython.display import SVG, display
except ImportError:
    warn('No Jupyter notebook installed. `.draw_nicely` will not work.')
    SVG = lambda *args, **kwargs: print('Install IPython...')
    display = lambda *args, **kwargs: print('Install IPython...')

try:
    import pymol2
except ImportError:
    warn('No Pymol module installed. `.make_pse` will not work.')
    pymol2 = None

########################################################################################################################


class _FragmensteinUtil:

    @classmethod
    def get_combined_rmsd(cls, followup_moved: Chem.Mol, followup_placed: Optional[Chem.Mol] = None,
                          hits: Optional[List[Chem.Mol]] = None) -> float:
        """
        The inbuilt RMSD calculations in RDKit align the two molecules, this does not align them.
        This deals with the case of multiple hits.
        For euclidean distance the square root of the sum of the differences in each coordinates is taken.
        For a regular RMSD the still-squared distance is averaged before taking the root.
        Here the average is done across all the atom pairs between each hit and the followup.
        Therefore, atoms in followup that derive in the blended molecule by multiple atom are scored multiple times.

        As a classmethod ``followup_placed`` and ``hits`` must be provided. But as an instance method they don't.

        :param followup_moved: followup compound moved by Igor or similar
        :param followup_placed: followup compound as placed by Fragmenstein
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
        mcs = rdFMCS.FindMCS([self.scaffold, self.initial_mol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrder)
        return Chem.MolFromSmarts(mcs.smartsString).GetNumAtoms()

    @property
    def percent_common(self) -> int:
        return round(self.num_common / self.initial_mol.GetNumAtoms() * 100)

    def stdev_from_mol(self, mol: Chem.Mol=None):
        """
        these values are stored from Fragmenstein for scaffold, chimera and positioned_mol

        :param mol: Chem.Mol
        :return: stdev list for each atom
        """
        if mol is None:
            mol = self.positioned_mol
        return [atom.GetDoubleProp('_Stdev') if atom.HasProp('_Stdev') else 0  for atom in mol.GetAtoms()]

    def origin_from_mol(self, mol: Chem.Mol = None):
        """
        these values are stored from Fragmenstein for scaffold, chimera and positioned_mol

        :param mol: Chem.Mol
        :return: stdev list for each atom
        """
        if mol is None:
            mol = self.positioned_mol
        if mol.HasProp('_Origins'):
            return json.loads(mol.GetProp('_Origins'))
        origin = []
        for atom in mol.GetAtoms():
            if atom.HasProp('_Origin'):
                x = atom.GetProp('_Origin')
                if x == 'none':
                    origin.append([])
                else:
                    origin.append(json.loads(x))
            else:
                origin.append([])
        return origin


    # class attribute for next method
    _i = 0

    def save_temp(self, mol):
        """
        This is a silly debug-by-print debug method. drop it in where you want to spy on stuff.
        """
        Chem.MolToMolFile(mol, f'debug_temp{self.i}.mol', kekulize=False)
        self._i += 1

    def save_commonality(self, filename:Optional[str]=None):
        """
        Saves an SVG of the followup fragmenstein with the common atoms with the chimeric scaffold highlighted.

        :param filename: optinal filename to save it as. Otherwise returns a Draw.MolDraw2DSVG object.
        :return:
        """
        mcs = rdFMCS.FindMCS([self.chimera, self.positioned_mol],
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

    def make_pse(self, filename='test.pse'):
        assert '.pse' in filename, 'Must be a pymol pse extension!'
        with pymol2.PyMOL() as pymol:
            tints = iter(['wheat', 'palegreen', 'lightblue', 'paleyellow', 'lightpink', 'palecyan', 'lightorange', 'bluewhite'])
            #pymol.cmd.bg_color('white')
            for h, hit in enumerate(self.hits):
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit, kekulize=False), f'hit{h}')
                pymol.cmd.color(next(tints), f'hit{h} and name C*')
            if self.scaffold:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.scaffold, kekulize=False), f'scaffold')
                pymol.cmd.color('tv_blue', f'scaffold and name C*')
            if self.chimera:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.chimera, kekulize=False), f'chimera')
                pymol.cmd.color('cyan', f'chimera and name C*')
            if self.positioned_mol:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.positioned_mol, kekulize=False), f'followup')
                pymol.cmd.color('tv_green', f'followup and name C*')
            if self.scaffold_options:
                for i, mol in enumerate(self.scaffold_options):
                    pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), f'opt{i}')
                    pymol.cmd.color('grey50', f'opt{i} and name C*')
            pymol.cmd.hide('sticks')
            pymol.cmd.hide('cartoon') # there should not be....
            pymol.cmd.show('lines', 'not polymer')
            if self.chimera:
                pymol.cmd.show('sticks', 'chimera')
            if self.positioned_mol:
                pymol.cmd.show('sticks', 'followup')
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
            #x = Chem.MolFromSmiles(Chem.MolToSmiles(x, kekuleSmiles=False), sanitize=False)
            Draw.PrepareAndDrawMolecule(d, x, **kwargs)
            d.FinishDrawing()
            if show:
                display(SVG(d.GetDrawingText()))
            return d
        except Exception as err:
            warn(f'*{err.__class__.__name__}* : {err}')
            display(x)



