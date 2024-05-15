import itertools
from typing import Union, Sequence, List, Iterator

import pandas as pd
import pebble
from rdkit import Chem, rdBase

from ._base import LabBench, binarize, unbinarize
from ..igor import pyrosetta  # this may be pyrosetta or a mock for Sphinx in RTD


class LabCombine(LabBench):

    def combine_subprocess(self, binary_hits: List[bytes]) -> dict:
        """
        This is the combination subprocess. The placement subprocess is ``place_subprocess``.
        They are very similar...
        """
        if self.Victor.uses_pyrosetta:
            pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        tentative_name = 'UNKNOWN'
        try:
            # this works if monster.fix_hits is not neeeded
            self.journal.debug(f'Combining {len(binary_hits)} hits')
            assert all([isinstance(hit, bytes) for hit in binary_hits]), 'Not all binary_hits are binary'
            hits: List[Chem.Mol] = [hit for hit in map(unbinarize, binary_hits) if hit is not None]
            assert len(hits) > 0, f'No valid hits ({len(binary_hits)} provided)'
            assert all([hit.GetNumAtoms() > 0 for hit in hits]), 'Some hits have no atoms!'
            if all([mol.HasProp('_Name') for mol in hits]):
                tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
                if tentative_name in self.blacklist:
                    raise ValueError(f'{tentative_name} is blacklisted')
            # `self.Victor` is likely `Victor` but the user may have switched for a subclass, cf. `VictorMock`...
            self.journal.debug(f'Using {self.Victor.__name__}')
            victor = self.Victor(hits=hits,
                                 pdb_block=self.pdbblock,
                                 ligand_resn='LIG',
                                 ligand_resi=self.ligand_resi,
                                 covalent_resi=self.covalent_resi,
                                 # a random residue is **still** required for the constaint ref atom.
                                 **self.settings
                                 )
            # the names may have been fixed by ``monster.fix_hits``
            tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
            if tentative_name in self.blacklist:
                raise ValueError(f'{tentative_name} is blacklisted')
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
                permute: bool = True,
                combination_size: int = 2,
                **kwargs) -> Union[pebble.ProcessMapFuture, pd.DataFrame]:
        """
        Combine all of ``mols`` with each other in combinations of ``combination_size``.
        Due to the way Monster works merging A with B may yield a different result to B with A.
        Hence the ``permute`` boolean argument.
        """  # extended at end of file.
        iterator: Iterator
        if permute:
            iterator = itertools.permutations(map(binarize, mols), combination_size)
        else:
            iterator = itertools.combinations(map(binarize, mols), combination_size)
        df = self(iterator=iterator, fun=self.combine_subprocess, **kwargs)
        df['outcome'] = df.apply(self.categorize, axis=1)
        with rdBase.BlockLogs():
            if 'unmin_binary' in df.columns:
                # doing at the binary level in case it failed
                df['simple_smiles'] = df.unmin_binary.apply(unbinarize)\
                                                     .apply(lambda m: m if m else Chem.Mol()) \
                                                     .apply(self.Victor.to_simple_smiles)
        self.fix_intxns(df)
        return df

    def twoway_combine(self,
                       primary_mols: Sequence[Chem.Mol],
                       secondary_mols: Sequence[Chem.Mol],
                       combination_size: int = 2,
                       **kwargs) -> Union[pebble.ProcessMapFuture, pd.DataFrame]:
        """
        Combine ``primary_mols`` with ``secondary_mols``.
        """
        iterator: Iterator
        if combination_size == 2:
            iterator = itertools.product(map(binarize, primary_mols), map(binarize, secondary_mols))
        elif combination_size > 2:
            extras = map(binarize, list(primary_mols) + list(secondary_mols))
            iterator = itertools.product(map(binarize, primary_mols), map(binarize, secondary_mols), extras,
                                         repeat=combination_size - 2)
        else:
            raise ValueError(f'combination_size must be > 2 (given: {combination_size}')
        df = self(iterator=iterator, fun=self.combine_subprocess, **kwargs)
        df['outcome'] = df.apply(self.categorize, axis=1)
        return df

    def serial_combine(self,
              hit_combinations: Sequence[Sequence[Chem.Mol]],
              **kwargs) -> Union[pebble.ProcessMapFuture, pd.DataFrame]:
        """
        This is akin to the way place works (table with hits)
        """  # extended at end of file.

        def generator():
            for hits in hit_combinations:
                yield [binarize(m) for m in hits]

        df = self(iterator=generator(), fun=self.combine_subprocess, **kwargs)
        df['outcome'] = df.apply(self.categorize, axis=1)
        return df


# prepend docstring to combine
LabCombine.combine.__doc__ = LabBench.__call__.__doc__ + '\n\n' + LabCombine.combine.__doc__
LabCombine.twoway_combine.__doc__ = LabBench.__call__.__doc__ + '\n\n' + LabCombine.twoway_combine.__doc__
