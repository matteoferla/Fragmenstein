from ._combine import LabCombine
from ._place import LabPlace

class Laboratory(LabCombine, LabPlace):
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
    """