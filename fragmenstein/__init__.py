########################################################################################################################
import os

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.5"
__citation__ = ""

########################################################################################################################

from warnings import warn


try:
    from .igor import Igor
except ImportError as err:
    warn(f'Igor (minimiser) unavailable — {err}.', category=ImportWarning)

try:
    from .victor import Victor
except ZeroDivisionError as err:
    warn(f'Victor (pipeline) unavailable — {err}.', category=ImportWarning)

from .monster import Monster
from .m_rmsd import mRSMD

FRAGMENSTEIN_DIR =  os.path.abspath(os.path.join(__file__, "../.."))
MOLS_EXAMPLES_DIR = os.path.join(FRAGMENSTEIN_DIR, "test_mols")