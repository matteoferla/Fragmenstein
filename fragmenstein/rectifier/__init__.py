from __future__ import annotations
########################################################################################################################

__doc__ = \
    """
Fix issue in auto-merging.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

from ._base import _RectifierBaseMixin  # provides the __init__ and shared methods
from ._ring import _RectifierRingMixin  # fixes rings
from ._odd import _RectifierOddMixin  # fixes specific oddities
from ._valence import _RectifierValenceMixin  # fixes valence
from rdkit import Chem


class Rectifier(_RectifierRingMixin, _RectifierOddMixin, _RectifierValenceMixin):
    """
    Fixes the nastiness.

    Do note that that the Chem.Mol does not get modified in place.
    ``.rwmol`` does and is a Chem.RWMol. The ``.mol`` is a Chem.Mol.

    The steps can be found in ``.modifications``.

    The .journal log is not given a handler.

    New atoms with have the bool prop ``_Novel``.

    Does not link distant atoms. For that see joining methods in Fragmenstein.

    >>> Rectifier(mol).fix().mol
    """

    def fix(self) -> Rectifier:
        self.fix_rings()  # from _RectifierRingMixin
        self.prevent_oddities()  # from _RectifierOddMixin
        self.ununspecified_bonds()  # from _RectifierValenceMixin
        self.triage_rings()  # from _RectifierValenceMixin
        Chem.Cleanup(self.rwmol)
        self.fix_issues()  # from _RectifierValenceMixin
        Chem.SanitizeMol(self.rwmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        return self
