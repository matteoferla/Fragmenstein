########################################################################################################################

__doc__ = \
    """
    ...
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################

from warnings import warn


try:
    from .egor import Egor
except ImportError as err:
    warn(f'Egor (minimiser) unavailable — {err}.', category=ImportWarning)

try:
    from .victor import Victor
except ZeroDivisionError as err:
    warn(f'Victor (pipeline) unavailable — {err}.', category=ImportWarning)

from .core import Fragmenstein
from .m_rmsd import mRSMD
