import logging
import pebble
import operator
from typing import (Any, Callable, Union, Iterator, Sequence, List)
from collections import Counter
import pandas as pd
from rdkit import Chem
from ..monster import Monster
from ..victor import Victor


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
    if isinstance(bin, float):  # nan
        return None
    if bin is None:
        return None
    try:
        return Chem.Mol(bin)
    except exception:
        return None

class LabBench:

    # the ``outcome`` column in the pandas dataframe can have these values in order of niceness:
    category_labels = ['crashed', 'too distant', 'timeout', 'unstable', 'equally sized', 'deviant', 'acceptable']

    def __init__(self, pdbblock: str, covalent_resi: Union[int, str, None]):
        self.pdbblock = pdbblock
        self.covalent_resi = covalent_resi
        self.init_options = '-ex1 -ex2 -no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false'
        self.raw_results = []
        if not len(Victor.journal.handlers):
            Victor.enable_stdout(logging.CRITICAL)

    # keep it just in case: (not called within class as that would require `self.__class__.binarize`)
    binarize: Callable[[Chem.Mol, bool], bytes] = staticmethod(binarize)
    unbinarize: Callable[[bytes, bool], Union[Chem.Mol, None]] = staticmethod(unbinarize)

    def get_completed(self, futures: pebble.ProcessMapFuture) -> pd.DataFrame:
        """
        A pebble.ProcessMapFuture cannot be stored as an instance attribute as it's not picklable
        as it handles the lock. As a result it must be passed.

        Fills ``self.raw_results`` with the results of the futures, before returning a dataframe.
        """
        result_iter: Iterator[int] = futures.result()
        self.raw_results = []
        while True:
            try:
                result = next(result_iter)
                self.raw_results.append(result)
            except TimeoutError as error:
                print("Function took longer than %d seconds" % error.args[1])
                self.raw_results.append({'error': 'TimeoutError', 'name': ''})
            except StopIteration as error:
                # it would be nice having a mock entry with the expected values...
                # self.raw_results.append({})
                break
            except KeyboardInterrupt as error:
                print('Keyboard!')
                self.raw_results.append({'error': 'KeyboardInterrupt', 'name': ''})
                break
            except Exception as error:
                print(f'{error.__class__.__name__}: {error}')
                self.raw_results.append({'error': error.__class__.__name__, 'name': ''})
        # list of dict to dataframe
        df = pd.DataFrame(self.raw_results)
        df['LE'] = df.apply(
            lambda row: row['∆∆G'] / (row.N_constrained_atoms + row.N_unconstrained_atoms),
            axis=1)
        nan_to_list = lambda value: value if isinstance(value, list) else []
        df['disregarded'] = df.disregarded.apply(nan_to_list)  # str
        df['regarded'] = df.regarded.apply(nan_to_list)  # str
        df['unminimized_mol'] = df.unmin_binary.apply(unbinarize)
        df['minimized_mol'] = df.min_binary.apply(unbinarize)
        df['hit_mols'] = df.hit_binaries.apply(lambda l: [unbinarize(b) for b in l] if isinstance(l, Sequence) else [])
        df['outcome'] = df.apply(self.categorize, axis=1)
        df['percent_hybrid'] = df.unminimized_mol.apply(self.percent_hybrid)
        return df

    def categorize(self, row: pd.Series) -> str:
        """
        Given a row categorise the 'outcome' field.
        Called by ``get_completed``.

        Subclass Laboratory to change what the outcome field ends up being.
        But do remember to update the ``category_labels`` attribute.
        """
        # see category_labels for list of values.
        is_filled: Callable[[Any], int] = lambda value: len(value) != 0 if hasattr(value, '__len__') else False
        if is_filled(row.disregarded) or 'ConnectionError' in row.error:
            # either is the same thing, but the Victor error catching makes them different
            return 'too distant'
        elif 'TimeoutError' in row.error:
            return 'timeout'
        elif is_filled(row.error):
            return 'crashed'
        elif max(map(operator.methodcaller('GetNumHeavyAtoms'), row.hit_mols))  \
                >= row.unminimized_mol.GetNumHeavyAtoms():
            return 'equally sized'
        elif row.comRMSD > 1:
            return 'too moved'
        elif row['∆∆G'] >= 0:
            return 'too contorted'
        else:
            return 'acceptable'

    def percent_hybrid(self, mol: Chem.Mol) -> float:
        """
        Given the origins how much of the molecule is not accounted for by the primary hit?

        :param mol:
        :return:
        """
        if not isinstance(mol, Chem.Mol):
            return float('nan')
        origins: List[List[str]] = Monster.origin_from_mol(None, mol)
        c = Counter([o[0].split('.')[0] if o else None for o in origins]).most_common()
        return round((mol.GetNumHeavyAtoms() - c[0][1]) / mol.GetNumHeavyAtoms() * 100, 1)

    def __call__(self,
                 iterator: Iterator,
                 fun: Callable,
                 n_cores: int = 4,
                 timeout: int = 240,
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
