########################################################################################################################

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2022 A.D."
__license__ = "MIT"
__version__ = "0.7"
__citation__ = ""

########################################################################################################################

from warnings import warn


try:
    from .igor import Igor
except ImportError as err:
    warn(f'Igor (minimizer) unavailable —{err}.', category=ImportWarning)

try:
    from .victor import Victor
except ImportError as err:
    warn(f'Victor (pipeline) unavailable —{err}.', category=ImportWarning)

try:
    from .multivictor import MultiVictorPlacement
except ImportError as err:
    warn(f'Victor (pipeline) unavailable —{err}.', category=ImportWarning)

from .monster import Monster
from .walton import Walton
from .laboratory import Laboratory
from molecular_rectifier import Rectifier
from .m_rmsd import mRMSD
from .mpro import MProVictor
from .mpro import data as mpro_data
from . import branding
