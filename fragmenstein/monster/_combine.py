########################################################################################################################

__doc__ = \
    """
Combine = merge/join
    """

########################################################################################################################

from rdkit import Chem
from rdkit.Chem import AllChem

from ._collapse_ring import _MonsterRing
from ._base import _MonsterBase
from ._communal import _MonsterCommunal
from ._merge import _MonsterMerge
from molecular_rectifier import Rectifier


########################################################################################################################

class _MonsterCombine(_MonsterRing, _MonsterMerge):

    def combine(self, keep_all: bool=True, collapse_rings: bool=True, joining_cutoff: int = 5):
        """
        Merge/links the hits. (Main entrypoint)

        :param keep_all:
        :param collapse_rings:
        :param joining_cutoff:
        :return:
        """
        # blanking
        self.unmatched = []
        self.modifications = {}
        # The following override class declared attributes.
        self.joining_cutoff = joining_cutoff
        self.throw_on_discard = keep_all
        # merge!
        if collapse_rings:
            col_hits = self.collapse_mols(self.hits)
            self.keep_copies(col_hits, 'Collapsed hits')
        else:
            col_hits = self.hits
        self.positioned_mol = self.simply_merge_hits(col_hits, linked=False)
        self.keep_copy(self.positioned_mol, 'merged template')
        ## Discard can happen for other reasons than disconnect
        if self.throw_on_discard and len(self.unmatched):
            raise ConnectionError(f'Could not combine with {self.unmatched} (>{self.joining_cutoff}')
        # expand and fix
        self.journal.debug(f'Merged')
        if collapse_rings:
            self.positioned_mol = self.expand_ring(self.positioned_mol)
        # bonded_as_original=False no longer needed.
        self.keep_copy(self.positioned_mol, 'expanded')
        self._join_internally(self.positioned_mol)
        self.journal.debug(f'Expanded')
        try:
            self.rectify()
        except RecursionError:
            self.journal.critical(f'Recursion limit in rectifier')
            raise ConnectionError(f'Can not rectify {self.positioned_mol}')
        self.journal.debug(f'Rectified')
        return self

    def rectify(self):
        recto = Rectifier(self.positioned_mol)
        try:
            recto.fix()
        except ConnectionError:
            self.journal.critical(f'This really odd cornercase: Rectifier broke the mol.')
            mol = self._emergency_joining(recto.mol)
            recto = Rectifier(mol)
            recto.fix()
        frags = AllChem.GetMolFrags(recto.mol, asMols=True)
        if len(frags) == 1:
            # all good
            self.positioned_mol = recto.mol
        else: # the molecule is still nasty. Getting largest.
            self.positioned_mol = sorted(frags, key=lambda mol: mol.GetNumAtoms(), reverse=True)[0]
        self.keep_copies(recto.modifications, 'fixed')
