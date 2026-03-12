__all__ = ['pyrosetta']

from importlib import util, machinery
from importlib.abc import MetaPathFinder
import sys
import warnings


class AttributeFilledMock:
    """
    This allows:

    .. code-block:: python
        mock = AttributeFilledMock()
        mock.me
        mock.me.again.over()

    It also fakes being a package so that
    ``import pyrosetta.rosetta.core`` doesn't explode.
    """
    warning_msg = 'PyRosetta is not installed! Yet it has been called. ' + \
                  'The mock object taking its place does nothing.' + \
                  'This is bound to raise an error'
    __signature__ = None
    __path__ = []
    __file__ = 'pyrosetta (mock)'
    __package__ = 'pyrosetta'

    def __init__(self, name: str = 'pyrosetta'):
        self.__name__ = name
        self.__spec__ = machinery.ModuleSpec(name, None, is_package=True)

    def __getattr__(self, attr: str):
        # ↓ create a child mock and register it in sys.modules
        child_name = f'{self.__name__}.{attr}'
        child = AttributeFilledMock(child_name)
        sys.modules[child_name] = child
        # ↓ cache on self so repeated access returns the same object
        object.__setattr__(self, attr, child)
        return child

    def __call__(self, *args, **kargs):
        warnings.warn('This call does nothing as PyRosetta is not installed',
                      category=RuntimeWarning)
        return self

    def __iter__(self):
        return iter([])

    def __repr__(self):
        return f'<AttributeFilledMock({self.__name__})>'


class _PyrosettaMockFinder(MetaPathFinder):
    """Intercepts any ``import pyrosetta.*`` and returns a mock,
    so ``import pyrosetta.rosetta.core as prc`` works without PyRosetta."""

    def find_module(self, fullname, path=None):
        if fullname == 'pyrosetta' or fullname.startswith('pyrosetta.'):
            return self
        return None

    def find_spec(self, fullname, path, target=None):
        if fullname == 'pyrosetta' or fullname.startswith('pyrosetta.'):
            return machinery.ModuleSpec(fullname, loader=self, is_package=True)
        return None

    def create_module(self, spec):
        return AttributeFilledMock(spec.name)

    def exec_module(self, module):
        sys.modules[module.__name__] = module


# ======================================================================================================================

if util.find_spec('pyrosetta'):
    import pyrosetta
else:
    warnings.warn('PyRosetta is not installed. A mock object is loaded. Any Igor calls will fail.',
                  category=RuntimeWarning)
    # ↓ install the finder *before* registering the root mock
    # so all future ``import pyrosetta.X.Y`` statements are intercepted
    sys.meta_path.insert(0, _PyrosettaMockFinder())
    pyrosetta = AttributeFilledMock()
    sys.modules['pyrosetta'] = pyrosetta
