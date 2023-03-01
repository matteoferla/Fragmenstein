from ._combine import LabCombine
from ._place import LabPlace
from ._base import binarize, unbinarize
from ._place import MolPlacementInput, BinPlacementInput
from ._extras import LabExtras
class Laboratory(LabCombine, LabPlace, LabExtras):
    """
    This class runs the combination or placement tasks of a list of molecules as subprocesses.
    The module used is ``pebble``, which is the same as ``multiprocessing`` but with a few more features.
    For more advanced usages, dask is a good candidate,
    but dask jobs on a cluster with a batch-queueing engine are non-trivial
    to set up.

    .. code-block:: python
        lab = Laboratory(pdbblock=template, covalent_resi=covalent_resi)
        results: pd.DataFrame = lab.combine(hits)
        px.histogram(results,
                     x='outcome',
                     category_orders={'outcome': lab.category_labels},
            title='Distribution of outcome')

    In the case of subclassed Victor, the user can create a `Laboratory` instance,
    and replace `.Victor` with their subclass of Victor.
    However, there are a few requirements, so please see ``combine_subprocess`` and ``place_subprocess``,
    where it gets called.
    Laboratory does not contain handy empty methods for the user to override like Victor (

    `post_igor_step`)
    """