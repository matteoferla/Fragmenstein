########################################################################################################################

__doc__ = \
    """
Combine = merge/join
    """

########################################################################################################################

from rdkit import Chem

from ._collapse_ring import _MonsterRing
from ._base import _MonsterBase
from ._communal import _MonsterCommunal
from ._merge import _MonsterMerge
from ..rectifier import Rectifier


########################################################################################################################

class _MonsterCombine(_MonsterRing, _MonsterMerge):

    def merge(self, keep_all=True, collapse_rings=True, joining_cutoff: int = 5):
        """
        Merge/links the hits. (Main entrypoint)

        :param keep_all:
        :param collapse_rings:
        :param joining_cutoff:
        :return:
        """
        # The following override class declared attributes.
        self.joining_cutoff = joining_cutoff
        self.throw_on_discard = keep_all
        # merge!
        if collapse_rings:
            col_hits = self.collapse_mols(self.hits)
            self.keep_copies(col_hits, 'Collapsed hits')
        else:
            col_hits = self.hits
        self.scaffold = self.simply_merge_hits(col_hits)
        self.keep_copy(self.scaffold, 'merged template')
        ## Discard can happen for other reasons than disconnect
        if keep_all and len(self.unmatched):
            raise ConnectionError(f'Could not combine with {self.unmatched} (>{self.joining_cutoff}')
        # expand and fix
        self.journal.debug(f'Merged')
        if collapse_rings:
            self.positioned_mol = self.expand_ring(self.scaffold)
        # bonded_as_original=False no longer needed.
        self.keep_copy(self.positioned_mol, 'expanded')
        self.journal.debug(f'Expanded')
        self.rectify()
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
        self.positioned_mol = recto.mol
        self.keep_copies(recto.modifications, 'fixed')
