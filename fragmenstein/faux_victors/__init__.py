"""
Victor alternatives not meant for production use, but for testing and debugging.

* Mictor: MCS hack for negative benchmarking: Victor with MCSMerger, i.e. no positional info used.
* Wictor: Victor without pyrosetta, i.e. no energy minimization.

To run in ``Laboratory`` class set the ``.Victor`` attribute to the desired ``Victor`` alternative.

... code-block:: python
    from framgenstein import Laboratory
    from framgenstein.faux_victors import Mictor  # MCS hack for negative benchmarking
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None)
    lab.Victor = Mictor
    combinations: pd.DataFrame = lab.combine(hits, n_cores=28)
    combinations.to_pickle(f'👾👾👾.p'))
"""

from .mcs_victor import Mictor, MCSMerger  #
from .no_pyrosetta import Wictor  # Victor without pyrosetta