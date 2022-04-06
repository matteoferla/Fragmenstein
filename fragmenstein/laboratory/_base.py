import itertools
import logging, pebble
import pandas as pd
from rdkit import Chem
from fragmenstein import Victor
from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD
import pyrosetta_help as ph
from typing import (Any, Callable, ClassVar, Final, ForwardRef, Generic, Literal, Optional, Protocol, Tuple, Type, TypeVar, Union, AbstractSet, ByteString, Container, ContextManager, Hashable, ItemsView, Iterable, Iterator, KeysView, Mapping, MappingView, MutableMapping, MutableSequence, MutableSet, Sequence, Sized, ValuesView, Awaitable, AsyncIterator, AsyncIterable, Coroutine, Collection, AsyncGenerator, AsyncContextManager, Reversible, SupportsAbs, SupportsBytes, SupportsComplex, SupportsFloat, SupportsIndex, SupportsInt, SupportsRound, ChainMap, Counter, Deque, Dict, DefaultDict, List, OrderedDict, Set, FrozenSet, NamedTuple, TypedDict, Generator, AnyStr, cast, final, get_args, get_origin, get_type_hints, NewType, no_type_check, no_type_check_decorator, NoReturn, overload, runtime_checkable, Text, TYPE_CHECKING)

def binarize(mol:Chem.Mol, ignore_errors:bool=True) -> bytes:
    # this is convoluted as None is a common error outcome with RDKit
    # so making it less cryptic
    if isinstance(mol, bytes):
        return mol   # messed up. it was bytes.
    elif not isinstance(mol, Chem.Mol): # likely None
        if ignore_errors:
            mol = Chem.Mol()    # noqa
        elif mol is None:
            raise TypeError('A Chem.Mol was expected by got None. '+
                            'RDKit returns None if it fails an operation.')
        else:
            raise TypeError('the molecule provide is not a Chem.Mol, '+
                            f'but a {type(mol).__name__}')
    exception = Exception if ignore_errors else ()
    try:
        return mol.ToBinary(propertyFlags=0b00010111)    # noqa it's an instance.
    except exception as error:
        return Chem.Mol().ToBinary(propertyFlags=0b00010111)  # noqa


def unbinarize(bin: bytes, ignore_errors:bool=True) -> Union[Chem.Mol, None]:
    exception = Exception if ignore_errors else ()
    if isinstance(bin, Chem.Mol):  # lol
        return bin
    try:
        return Chem.Mol(bin)
    except exception:
        return None

class LabBench:

    # the ``outcome`` column in the pandas dataframe can have these values in order of niceness:
    category_labels = ['crashed', 'too distant', 'timeout', 'unstable', 'deviant', 'acceptable']

    def __init__(self, pdbblock: str, covalent_resi: Union[int, str, None]):
        self.pdbblock = pdbblock
        self.covalent_resi = covalent_resi
        self.init_options = '-no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false'

    # keep it just in case: (not called within class as that would require `self.__class__.binarize`)
    binarize: Callable[[Chem.Mol, bool], bytes] = staticmethod(binarize)
    unbinarize: Callable[[bytes, bool], Union[Chem.Mol, None]] = staticmethod(unbinarize)

    def get_completed(self, futures: pebble.ProcessMapFuture) -> pd.DataFrame:
        """
        A pebble.ProcessMapFuture cannot be stored as an instance attribute as it's not picklable
        as it handles the lock. As a result it must be passed.
        """
        result_iter: Iterator[int] = futures.result()
        results = []
        while True:
            try:
                result = next(result_iter)
                results.append(result)
            except TimeoutError as error:
                print("Function took longer than %d seconds" % error.args[1])
                results.append({'error': 'TimeoutError', 'name': ''})
            except StopIteration as error:
                results.append({})
                break
            except KeyboardInterrupt as error:
                print('Keyboard!')
                results.append({'error': 'KeyboardInterrupt', 'name': ''})
                break
            except Exception as error:
                print(f'{error.__class__.__name__}: {error}')
                results.append({'error': error.__class__.__name__, 'name': ''})
        # list of dict to dataframe
        combinations = pd.DataFrame(results)
        combinations['LE'] = combinations.apply(
            lambda row: row['∆∆G'] / (row.N_constrained_atoms + row.N_unconstrained_atoms),
            axis=1)
        nan_to_list = lambda value: value if isinstance(value, list) else []
        combinations['disregarded'] = combinations.disregarded.apply(nan_to_list)
        combinations['regarded'] = combinations.regarded.apply(nan_to_list)
        combinations['unminimized_mol'] = combinations.unmin_binary.apply(unbinarize)
        combinations['minimized_mol'] = combinations.min_binary.apply(unbinarize)

        def categorize(row: pd.Series) -> str:
            # see category_labels for list of values.
            is_filled: Callable[[Any], int] = lambda value: len(value) != 0 if hasattr(value, '__len__') else False
            if is_filled(row.disregarded):
                return 'too distant'
            elif row.error == 'TimeoutError':
                return 'timeout'
            elif is_filled(row.error):
                return 'crashed'
            elif row.comRMSD > 1:
                return 'too moved'
            elif row['∆∆G'] >= 0:
                return 'too contorted'
            else:
                return 'acceptable'

        combinations['outcome'] = combinations.apply(categorize, axis=1)
        return combinations

    def __call__(self,
                 iterator: Iterator,
                 fun: Callable,
                 n_cores: int = 4,
                 timeout: int = 120,
                 asynchronous: bool = False
                 ):
        """Combine/permute the molecules ``mols``
        on ``n_cores`` subprocesses.
        killing any that live longer than ``timeout`` seconds.
        The method returns an iterator of promises ``pebble.ProcessMapFuture`` if ``asynchronous`` is True,
        or the results as a pandas DataFrame. To convert the promises to a dataframe use ``get_completed``."""
        with pebble.ProcessPool(max_workers=n_cores, max_tasks=n_cores) as pool:
            futures: pebble.ProcessMapFuture = pool.map(fun, iterator, timeout=timeout)
        if asynchronous:
            return futures
        else:
            return self.get_completed(futures)
