import pebble
from typing import (Any, Callable, Union, Iterator)

import pandas as pd
import pebble
from rdkit import Chem


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
        df = pd.DataFrame(results)
        df['LE'] = df.apply(
            lambda row: row['∆∆G'] / (row.N_constrained_atoms + row.N_unconstrained_atoms),
            axis=1)
        nan_to_list = lambda value: value if isinstance(value, list) else []
        df['disregarded'] = df.disregarded.apply(nan_to_list)
        df['regarded'] = df.regarded.apply(nan_to_list)
        df['unminimized_mol'] = df.unmin_binary.apply(unbinarize)
        df['minimized_mol'] = df.min_binary.apply(unbinarize)
        df['hit_mols'] = df.hit_binaries.apply(lambda l: [unbinarize(b) for b in l])


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
