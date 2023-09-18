import logging
import pebble
import operator
from typing import (Any, Callable, Union, Iterator, Sequence, List, Dict)
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from ..monster import Monster
from ..victor import Victor

# not needed for binarize... but just in case user is not using them...
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

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
    except KeyboardInterrupt as err:
        raise err
    except KeyboardInterrupt as err:
        raise err
    except exception as error:
        return Chem.Mol().ToBinary(propertyFlags=0b00010111)  # noqa


def unbinarize(bin: bytes, ignore_errors:bool=True) -> Union[Chem.Mol, None]:
    exception = Exception if ignore_errors else ()
    if isinstance(bin, Chem.Mol):  # lol
        return bin
    elif not isinstance(bin, bytes) and ignore_errors:
        return None
    elif not isinstance(bin, bytes):
        raise TypeError('the molecule binary provide is not a bytes')
    try:
        return Chem.Mol(bin)
    except KeyboardInterrupt as err:
        raise err
    except exception:
        return None

class LabBench:

    # the ``outcome`` column in the pandas dataframe can have these values in order of niceness:
    category_labels = ['crashed', 'too distant', 'timeout', 'unstable', 'equally sized', 'deviant', 'acceptable']
    Victor = Victor  # So it can be swapped for a subclass w/o the need to subclass Laboratory

    def __init__(self, pdbblock: str,
                 covalent_resi: Union[int, str, None] = None,
                 ligand_resi: Union[str, None]=None,
                 run_plip: bool=False,
                 ):
        self.pdbblock = pdbblock
        self.covalent_resi = covalent_resi
        self.init_options = '-ex1 -ex2 -no_optH false -mute all -ignore_unrecognized_res true -load_PDB_components false'
        self.raw_results = []
        self.ligand_resi = ligand_resi
        self.run_plip = run_plip
        self.blacklist = []  # list of names to skip
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
                Victor.journal.error("Function took longer than %d seconds" % error.args[1])
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
                Victor.journal.error(f'{error.__class__.__name__}: {error}')
                self.raw_results.append({'error': error.__class__.__name__, 'name': ''})
        # list of dict to dataframe
        df = pd.DataFrame(self.raw_results)
        if not len(df) or '∆∆G' not in df.columns:
            Victor.journal.critical('No results were found. Returning an empty dataframe.')
            return df
        df['LE'] = df.apply(
            lambda row: - row['∆∆G'] / (row.N_constrained_atoms + row.N_unconstrained_atoms),
            axis=1)
        nan_to_list = lambda value: value if isinstance(value, list) else []
        nan_to_mol = lambda m: m if isinstance(m, Chem.Mol) else Chem.Mol()
        df['disregarded'] = df.disregarded.apply(nan_to_list)  # str
        df['regarded'] = df.regarded.apply(nan_to_list)  # str
        df['unminimized_mol'] = df.unmin_binary.apply(unbinarize).apply(nan_to_mol)
        df['minimized_mol'] = df.min_binary.apply(unbinarize).apply(nan_to_mol)
        df['hit_mols'] = df.hit_binaries.apply(lambda l: [unbinarize(b) for b in l] if isinstance(l, Sequence) else [])
        df['hit_names'] = df.hit_mols.apply(lambda v: [m.GetProp('_Name') for m in v])
        df['percent_hybrid'] = df.unminimized_mol.apply(self.percent_hybrid)
        # if plipped fix nan interactions
        self.fix_intxns(df)
        # macrocyclics... yuck.
        df['largest_ring'] = df.minimized_mol.apply(lambda mol: max([0] + list(map(len, mol.GetRingInfo().AtomRings()))))
        df['N_HA'] = df.minimized_mol.apply(Chem.Mol.GetNumHeavyAtoms)
        df['N_rotatable_bonds'] = df.minimized_mol.apply(AllChem.CalcNumRotatableBonds)
        return df

    def categorize(self,
                   row: pd.Series,
                   size_tolerance: int=0,
                   move_cutoff:float=1.,
                   ddG_cutoff:float=0.,
                   ) -> str:
        """
        Given a row categorise the 'outcome' field.
        Called by ``get_completed``.

        size_tolerance is the number of atoms below the largest hit that is still acceptable.
        move_cutoff is the RMSD of the minimized molecules. Below this is acceptable.
        ddG_cutoff is the ∆∆G of the minimized molecules. Below this is acceptable.

        Subclass Laboratory to change what the outcome field ends up being.
        But do remember to update the ``category_labels`` attribute.
        """
        # see category_labels for list of values.
        is_filled: Callable[[Any], int] = lambda value: len(value) != 0 if hasattr(value, '__len__') else False
        if 'disregarded' not in row.index:
            return 'crashed'  # the whole row is empty
        elif is_filled(row.disregarded) or 'DistanceError' in row.error:
            # either is the same thing, but the Victor error catching makes them different
            return 'too distant'
        elif 'TimeoutError' in row.error:
            return 'timeout'
        elif is_filled(row.error):
            return 'crashed'
        elif max(map(operator.methodcaller('GetNumHeavyAtoms'), row.hit_mols)) - size_tolerance  \
                >= row.unminimized_mol.GetNumHeavyAtoms():
            return 'equally sized'
        elif row.comRMSD > move_cutoff:
            return 'too moved'
        elif row['∆∆G'] >= ddG_cutoff:
            return 'too contorted'
        else:
            return 'acceptable'

    def percent_hybrid(self, mol: Chem.Mol) -> float:
        """
        Given the origins how much of the molecule is solely from the second hit?
        It's the ratio of the number of atoms from the second hit to the total number of atoms with an inspiration.
        Do note that monster has the dynamic attributes, ``.percent_common`` and ``.num_common``.
        This is how-much and how-many of the molecule is common to both hits.
        This metric is not in Monster but in the Laboratory, making it handy for overriding.
        It also aims at not diagnosing how much overlap there is but to diagnose how combined the molecule is.

        Basically there are a bunch of Qs to ask if there are two hits:

        * How much of a molecule is common to both?
        * How much is from either single hit?
        * How much is unmapped?

        But when discussing mergers,
        if a bunch of atoms are common, the first contribute more,
        but the second hit does not contribute anything
        then it's not a merger: it's just the first hit.
        Hence the need for this odd ratio.
        """
        if not isinstance(mol, Chem.Mol):
            return float('nan')
        # origins: list of one entry for each followup atom to list of strings of hit name dot atom name
        origins: List[List[str]] = Monster.origin_from_mol(None, mol)
        hit_names = [[ori_name.split('.')[0] for ori_name in atomic if ori_name] if atomic else [] for atomic in origins]
        flat_names = [name for atomic in hit_names for name in atomic]
        names = set(flat_names)
        if len(names) == 0:  # no inspirations
            return float('nan')
        elif len(names) == 1:  # one inspirations
            return 0.
        # hit to number of atoms with only that hit as origin
        single_origin: Dict[str, int] = {name: sum([name in atomic for atomic in hit_names if len(atomic) == len(names) - 1]) for name in names}
        sorted_names = sorted(single_origin, key=single_origin.get, reverse=True)
        if sum(single_origin.values()) == 0:
            return 0.
        return 100 - int(single_origin[sorted_names[0]] / sum(single_origin.values()) * 100)

    def __call__(self,
                 iterator: Iterator,
                 fun: Callable,
                 n_cores: int = 4,
                 timeout: int = 240,
                 max_tasks: int = 0,  # 0 mean infinity
                 asynchronous: bool = False
                 ):
        """Combine/permute the molecules ``mols``
        on ``n_cores`` subprocesses.
        killing any that live longer than ``timeout`` seconds.
        The method returns an iterator of promises ``pebble.ProcessMapFuture`` if ``asynchronous`` is True,
        or the results as a pandas DataFrame. To convert the promises to a dataframe use ``get_completed``."""

        def max_out(inner_iterator, maximum: int):
            for i, item in zip(range(maximum), inner_iterator):
                yield item

        if max_tasks > 0:
            iterator = max_out(iterator, max_tasks)

        with pebble.ProcessPool(max_workers=n_cores, max_tasks=n_cores) as pool:
            futures: pebble.ProcessMapFuture = pool.map(fun, iterator, timeout=timeout)
        if asynchronous:
            return futures
        else:
            return self.get_completed(futures)

    @staticmethod
    def fix_intxns(df):
        """
        The interactions from PLIP will be columns with a list of tuples,
        of the form ('hbond', 'LEU', 127),
        _i.e._ interaction name, residue name and residue number.
        """
        intxn_names = [c for c in df.columns if isinstance(c, tuple)]
        for intxn_name in intxn_names:
            df[intxn_name] = df[intxn_name].fillna(0).astype(int)