import sys
from typing import (Union, Iterator, Sequence, List)

import pandas as pd
import pebble
from rdkit import Chem


from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD

from typing import TypedDict  # monkeypatched by .legacy. Absent in <Py3.8


from ._base import LabBench, binarize, unbinarize

class MolPlacementInput(TypedDict):
    """
    The place command accepts a list of dictionary or a DataFrame,
    with three keys.

    The problem is that hits needs to be pickled properly without losing Properties.
    Hence the conversion to BinPlacementInput

    """
    smiles: str
    name: str
    hits: Sequence[Chem.Mol]

class BinPlacementInput(TypedDict):
    """
    Like MolPlacementInput but with binary hits

    """
    smiles: str
    name: str
    binary_hits: Sequence[bytes]

class LabPlace(LabBench):

    def place_subprocess(self, inputs: BinPlacementInput):
        """
        This is the placement subprocess. The combination subprocess is ``place_subprocess``.
        They are very similar...
        """
        name: str = inputs['name']
        smiles: str = inputs['smiles']
        pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        try:
            binary_hits = inputs['binary_hits']
            hits: List[Chem.Mol] = [hit for hit in map(unbinarize, binary_hits) if hit]
            assert len(hits) > 0, 'No valid hits!'
            # `self.Victor` is likely `Victor` but the user may have switched for a subclass, cf. `VictorMock`...
            v = self.Victor(hits=hits,
                            pdb_block=self.pdbblock,
                            ligand_resn='LIG',
                            ligand_resi='1B',
                            covalent_resi=self.covalent_resi,
                            )
            v.place(smiles, long_name=name)
            result: dict = v.summarize()
            result['unmin_binary'] = binarize(v.monster.positioned_mol)
            result['min_binary'] = binarize(v.minimized_mol)
            result['hit_binaries'] = [binarize(h) for h in v.hits]
            return result
        except KeyboardInterrupt as err:
            raise err
        except Exception as error:
            error_msg = f'{error.__class__.__name__} {error}'
            v.journal.critical(f'*** {error_msg} for {name}')
            return dict(error=error_msg, name=name)

    def place(self,
              queries: Union[pd.DataFrame, Sequence[MolPlacementInput]],
              expand_isomers:bool=False,
              **kwargs) -> Union[pebble.ProcessMapFuture, pd.DataFrame]:
        """
        Due to the way Monster works merging A with B may yield a different result to B with A.
        Hence the ``permute`` boolean argument.
        """  # extended at end of file.
        if isinstance(queries, pd.DataFrame):
            assert 'smiles' in queries.columns
            assert 'name' in queries.columns
            assert 'hits' in queries.columns
            pre_iterator = queries.iterrows()
        elif isinstance(queries, Sequence):
            pre_iterator = enumerate(queries)
        else:   # it may crash...
            pre_iterator = queries
        iterator: Iterator[BinPlacementInput] = iter([{'smiles': d['smiles'],
                                                      'name': d['name'],
                                                      'binary_hits': [binarize(m) for m in d['hits']]}
                                                     for _, d in pre_iterator])
        if expand_isomers:
            iterator = [{'binary_hits': data['binary_hits'],
                         'name': data['name']+f'-isomer_{i}',
                         'smiles': sub_smiles}
                                for data in iterator
                                for i, sub_smiles in enumerate(self.Victor.get_isomers_smiles(data['smiles']))
                        ]
        return self(iterator=iterator, fun=self.place_subprocess, **kwargs)

# prepend docstring to combine
LabPlace.place.__doc__ = LabBench.__call__.__doc__ + '\n\n' + LabPlace.place.__doc__
