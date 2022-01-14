__all__ = ['pyrosetta']

from importlib import util
import sys
from warnings import warn

class AttributeFilledMock:
    """
    This allows:

    .. code-block:: python
        mock = AttributeFilledMock()
        mock.me
        mock.me.again.over()
    """

    warning_msg = 'PyRosetta is not installed! Yet it has been called. '+\
                  'The mock object taking its place does nothing.'+\
                  'This is bound to raise an error'

    def __getattr__(self, attr: str):
        return self

    def __call__(self, *args, **kargs):
        warn()
        return self

# ======================================================================================================================

if util.find_spec('pyrosetta'):
    import pyrosetta
else:
    warn('PyRosetta is not installed. A mock object is loaded. Any calls will fail.')
    pyrosetta = AttributeFilledMock()
    sys.modules['pyrosetta'] = pyrosetta

