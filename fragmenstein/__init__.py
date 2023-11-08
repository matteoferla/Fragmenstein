########################################################################################################################

__doc__ = \
    """
    See GitHub documentation
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2022 A.D."
__license__ = "MIT"
__citation__ = ""
from .version import __version__

########################################################################################################################

from warnings import warn
from . import legacy as _  # monkeypatches singledispatchmethod TypedDict Unpack on older Pythons

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

try:
    from .laboratory import Laboratory, binarize, unbinarize, place_input_validator
except ImportError as err:
    warn(f'Laboratory unavailable —{err}.', category=ImportWarning)

try:
    from .openmm import OpenVictor, Fritz
except ImportError as err:
    warn(f'OpenVictor / Fritz unavailable —{err}.', category=ImportWarning)

from .monster import Monster, MinizationOutcome
from .walton import Walton
from molecular_rectifier import Rectifier
from .m_rmsd import mRMSD
from .mpro import MProVictor
from .mpro import data as mpro_data
from . import branding
from .display import display_mols, MolNGLWidget
from .error import FragmensteinError, DistanceError, ShoddyCodeError, PoisonError
from .faux_victors import Wictor, Quicktor  # the rest are mostly experiments
from .settings import default_settings, cli_default_settings

if __name__ == '__main__':
    from .cli import main
    main()
