import itertools
from typing import (Union, Sequence, List)

import pandas as pd
import pebble
from rdkit import Chem

from ._base import LabBench, binarize, unbinarize
from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD


class LabCombine(LabBench):

    def combine_subprocess(self, binary_hits: List[bytes]):
        """
        This is the combination subprocess. The placement subprocess is ``place_subprocess``.
        They are very similar...
        """
        pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        tentative_name = 'UNKNOWN'
        try:
            hits: List[Chem.Mol] = [hit for hit in map(unbinarize, binary_hits) if hit]
            assert len(hits) > 0, 'No valid hits!'
            tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
            # `self.Victor` is likely `Victor` but the user may have switched for a subclass, cf. `VictorMock`...
            victor = self.Victor(hits=hits,
                       pdb_block=self.pdbblock,
                       ligand_resn='LIG',
                       ligand_resi=self.ligand_resi,
                       covalent_resi=self.covalent_resi,
                       # a random residue is **still** required for the constaint ref atom.
                       )
            victor.monster_throw_on_discard = True
            victor.monster.throw_on_discard = True
            victor.combine()
            result: dict = victor.summarize()
            result['unmin_binary'] = binarize(victor.monster.positioned_mol)
            result['min_binary'] = binarize(victor.minimized_mol)
            result['hit_binaries'] = [binarize(h) for h in victor.hits]
            if self.run_plip:
                result.update(victor.get_plip_interactions())
            return result
            # victor.make_pse()
        except KeyboardInterrupt as err:
            raise err
        except Exception as error:
            error_msg = f'{error.__class__.__name__} {error}'
            self.Victor.journal.critical(f'*** {error_msg} for {tentative_name}')
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
        df = self(iterator=iterator, fun=self.combine_subprocess, **kwargs)
        df['outcome'] = df.apply(self.categorize, axis=1)
        return df

# prepend docstring to combine
LabCombine.combine.__doc__ = LabBench.__call__.__doc__  + '\n\n' + LabCombine.combine.__doc__
