import itertools
from typing import (Union, Sequence, List)

import pandas as pd
import pebble
from rdkit import Chem

from fragmenstein import Victor
from ._base import LabBench, binarize
from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD


class LabCombine(LabBench):

    def combine_subprocess(self, binary_hits: List[bytes]):
        pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        tentative_name = 'UNKNOWN'
        try:
            hits: List[Chem.Mol] = [Chem.Mol(bh) for bh in binary_hits]
            tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
            v = Victor(hits=hits,
                       pdb_block=self.pdbblock,
                       ligand_resn='LIG',
                       ligand_resi='1B',
                       covalent_resi=self.covalent_resi,
                       # a random residue is **still** required for the constaint ref atom.
                       )
            v.combine()
            result: dict = v.summarize()
            result['unmin_binary'] = binarize(v.monster.positioned_mol)
            result['min_binary'] = binarize(v.minimized_mol)
            result['hit_binaries'] = [binarize(h) for h in v.hits]
            return result
            # v.make_pse()
        except Exception as error:
            error_msg = f'{error.__class__.__name__} {error}'
            Victor.journal.critical(f'*** {error_msg} for {tentative_name}')
            return dict(error=error_msg, name=tentative_name)

    def combine(self,
                mols: Sequence[Chem.Mol],
                permute:bool=True,
                combination_size:int=2,
                **kwargs) -> Union[pebble.ProcessMapFuture, pd.DataFrame]:
        """
        Due to the way Monster works merging A with B may yield a different result to B with A.
        Hence the ``permute`` boolean argument.
        """  # extended at end of file.
        if permute:
            iterator = itertools.permutations(map(binarize, mols), combination_size)
        else:
            iterator = itertools.combinations(map(binarize, mols), combination_size)
        return self(iterator=iterator, fun=self.combine_subprocess, **kwargs)

# prepend docstring to combine
LabCombine.combine.__doc__ = LabBench.__call__.__doc__  + '\n\n' + LabCombine.combine.__doc__
