import itertools
import logging, pebble
import pandas as pd
from rdkit import Chem
from fragmenstein import Victor
from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD
import pyrosetta_help as ph
from typing import (Any, Callable, ClassVar, Final, ForwardRef, Generic, Literal, Optional, Protocol, Tuple, Type, TypeVar, Union, AbstractSet, ByteString, Container, ContextManager, Hashable, ItemsView, Iterable, Iterator, KeysView, Mapping, MappingView, MutableMapping, MutableSequence, MutableSet, Sequence, Sized, ValuesView, Awaitable, AsyncIterator, AsyncIterable, Coroutine, Collection, AsyncGenerator, AsyncContextManager, Reversible, SupportsAbs, SupportsBytes, SupportsComplex, SupportsFloat, SupportsIndex, SupportsInt, SupportsRound, ChainMap, Counter, Deque, Dict, DefaultDict, List, OrderedDict, Set, FrozenSet, NamedTuple, TypedDict, Generator, AnyStr, cast, final, get_args, get_origin, get_type_hints, NewType, no_type_check, no_type_check_decorator, NoReturn, overload, runtime_checkable, Text, TYPE_CHECKING)


from ._base import LabBench, binarize, unbinarize

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
