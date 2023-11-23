import pandas as pd
from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import PandasTools
from .validator import place_input_validator
from .._cli_defaults import cli_default_settings
import time
from typing import List, Dict, Any, Optional


class LabExtras:
    error_classifications = {'SUCCESS': 'success',
                             'UNCLASSIFIED': 'unclassified',
                             'ConnectionError': 'distance',  # legacy
                             'DistanceError': 'distance',
                             'RectificationError': 'incorrect rectification',
                             'FullOverlapError': 'distance',
                             'TimeoutError': 'timeout',
                             'No sub-structure match': 'incorrect parameterisation #1',
                             'Violation occurred': 'incorrect rectification #1',
                             'KekulizeException': 'incorrect rectification #1',
                             'AtomValenceException': 'incorrect rectification #1',
                             'Bad Conformer Id': 'incorrect rectification #1',
                             ('rosetta', 'not found'): 'incorrect parameterisation #2',
                             ('rosetta', 'nan'): 'incorrect parameterisation #3',
                             'but not the minimised': 'incorrect parameterisation #4',
                             'UtilityExitException': 'incorrect parameterisation #5',
                             'TypeError': 'embedding error',
                             'RecursionError': 'embedding error',
                             'utterly impossible': 'embedding error'
                             }

    @classmethod
    def error_classify(cls, text):
        if not text:
            return cls.error_classifications['SUCCESS']
        text = str(text)
        for key, value in cls.error_classifications.items():
            if isinstance(key, tuple):
                if all(x in text for x in key):
                    return value
            elif key in text:
                return value
            else:
                return cls.error_classifications['UNCLASSIFIED']

    @classmethod
    def convert_to_sdf(self,
                       df: pd.DataFrame,
                       filename: str = f'fragmenstein.sdf',
                       acceptable_only=True,
                       sort_values: str = cli_default_settings['ranking'],
                       name_col: str = 'name',
                       mol_col='minimized_mol'):
        if acceptable_only:
            df = df.loc[(df.outcome == 'acceptable')]
        short = df.sort_values(sort_values).reset_index().drop_duplicates(name_col)

        PandasTools.WriteSDF(df=short,
                             out=filename,
                             molColName=mol_col,
                             idName=name_col,
                             properties=['regarded', 'smiles', '∆∆G', '∆G_bound', '∆G_unbound',
                                         'comRMSD', 'N_constrained_atoms', 'N_unconstrained_atoms', 'runtime',
                                         'LE', 'outcome',
                                         'percent_hybrid']
                             )

    # ---------- CLI --------------------------------------------

    @classmethod
    def core_ops(cls, hit_replacements, sw_databases, **settings):
        combinations: pd.DataFrame = cls._combine_ops(**settings)
        cls.correct_weaklings(hit_replacements, combinations)
        uncat_analogs: List[pd.DataFrame] = []
        for sw_db in sw_databases:
            s = cls.sw_search(combinations, sw_db=sw_db, **settings)
            if len(s):
                uncat_analogs.append(s)
        assert uncat_analogs, 'No analogues were found!'
        analogs: pd.DataFrame = pd.concat(uncat_analogs, ignore_index=True).drop_duplicates('smiles')
        analogs.to_pickle(f'fragmenstein_analogs{settings["suffix"]}.pkl.gz')
        placements: pd.DataFrame = cls._place_ops(analogs=analogs, **settings)
        return placements

    @classmethod
    def _combine_ops(cls,
                     hits,
                     pdbblock,
                     suffix: str = cli_default_settings['suffix'],
                     n_cores: int = cli_default_settings['n_cores'],
                     combination_size: int = cli_default_settings['combination_size'],
                     timeout: int = cli_default_settings['timeout'],
                     max_tasks: int = cli_default_settings['max_tasks'],
                     blacklist: List[str] = cli_default_settings['blacklist'],
                     **settings) -> pd.DataFrame:
        """
        One of the operations of ``core_ops``.
        A thin wrapper around ``combine``.
        """
        lab = cls(pdbblock=pdbblock, covalent_resi=None)  # noqa it's inherited later
        lab.blacklist = blacklist
        tick = time.time()
        combinations: pd.DataFrame = lab.combine(hits,  # noqa it's inherited later
                                                 n_cores=n_cores,
                                                 timeout=timeout,
                                                 combination_size=combination_size,
                                                 max_tasks=max_tasks)
        combinations.to_pickle(f'fragmenstein_mergers{suffix}.pkl.gz')
        combinations.to_csv(f'fragmenstein_mergers{suffix}.csv')
        print(tick - time.time())
        return combinations

    @classmethod
    def sw_search(cls, combinations: pd.DataFrame, suffix: str,
                  sw_dist: int, sw_length: int, top_mergers: int,
                  ranking: str, sw_db: str, ranking_ascending: Optional[bool] = None,
                  sws: Optional = None,
                  **setting) -> pd.DataFrame:
        if ranking_ascending is None:
            ranking_ascending = False if ranking in ('LE', 'N_interactions') else True
        queries = combinations.sort_values(ranking, ascending=ranking_ascending) \
            .loc[(combinations.outcome == 'acceptable')] \
            .drop_duplicates('simple_smiles') \
            .reset_index() \
            .head(top_mergers)
        if sws is None:
            from smallworld_api import SmallWorld
            sws = SmallWorld()
        analogs = sws.search_many(queries.simple_smiles.to_list(),
                                  dist=sw_dist,
                                  length=sw_length,
                                  db=sw_db,
                                  tolerated_exceptions=Exception)
        print(f'Found {len(analogs)} analogues')
        # query_index was added clientside to keep track!
        analogs['catalogue'] = sw_db
        analogs['query_name'] = analogs.query_index.map(queries.name.to_dict())
        analogs['hits'] = analogs.query_index.map(queries.hit_mols.to_dict())
        analogs['hit_names'] = analogs.hits.apply(lambda m: [mm.GetProp('_Name') for mm in m])
        analogs['minimized_merger'] = analogs.query_index.map(queries.minimized_mol.to_dict())
        analogs['unminimized_merger'] = analogs.query_index.map(queries.unminimized_mol.to_dict())
        analogs['name'] = analogs['id'] + ':' + analogs['query_name']
        analogs['smiles'] = analogs.hitSmiles.str.split(expand=True)[0]
        analogs['custom_map'] = analogs.apply(cls.get_custom_map, axis=1)
        analogs.to_pickle(f'fragmenstein_analogues{suffix}.{sw_db}.pkl.gz')
        return analogs

    @classmethod
    def _place_ops(cls, analogs, pdbblock, n_cores, timeout, suffix, **settings) -> pd.DataFrame:
        """
        This is the classmethod called by ``core_ops``.
        The instance method ``place`` does the actual work, this is a thin wrapper.
        """
        lab = cls(pdbblock=pdbblock, covalent_resi=None, run_plip=True)  # noqa it's inherited later
        placements: pd.DataFrame = lab.place(place_input_validator(analogs),
                                             n_cores=n_cores,
                                             timeout=timeout)
        placements.to_pickle(f'fragmenstein_placed{suffix}.pkl.gz')
        placements.to_csv(f'fragmenstein_placed{suffix}.csv')
        # print(placements.outcome.value_counts())
        return placements

    @classmethod
    def get_custom_map(cls, row: pd.Series) -> Dict[str, Dict[int, int]]:
        """
        SmallWorld returns a mapping of the indices.
        This functions maps the indices back to the original hit molecules.

        :param row: SW pd.Series (or dict)
        :return:
        """
        # faux PDB block to trick the safeguards against bad PDB blocks/filenames
        temp = cls.Victor(row.hits, pdb_block=Chem.MolToPDBBlock(Chem.MolFromFASTA('A')))
        temp.monster.positioned_mol = row.unminimized_merger
        temp.minimized_mol = row.minimized_merger
        return temp.migrate_sw_origins(row)

    @classmethod
    def correct_weaklings(cls, hit_replacements: pd.DataFrame, target_df: pd.DataFrame):
        """
        Some hits are weaker than the original fragment.
        If compounds in ``target_df`` are weaker than the inspiration in ``hit_replacements``,
        then the outcome is changed to ``weaker``.
        """
        # just in case...
        target_df['hit_names'] = target_df.hit_mols \
            .apply(lambda v: v if isinstance(v, list) else []) \
            .apply(lambda v: [m.GetProp('_Name') for m in v])
        # assign by name
        dG_mapping = hit_replacements.set_index('name')['∆∆G'].to_dict()
        dG_mapping.update({k.replace('-', '_'): v for k, v in dG_mapping.items()})
        get_lowest = lambda names: min([0] + [dG_mapping.get(name, 0) for name in names])
        target_df['lowest_hit_∆∆G'] = target_df.hit_names.apply(get_lowest)
        worseness_mask = (target_df['∆∆G'] > target_df['lowest_hit_∆∆G'] * 0.8) & (
                target_df.outcome == 'acceptable')
        target_df.loc[worseness_mask, 'outcome'] = 'weaker'

    @classmethod
    def replace_hits(cls,
                     pdbblock:str,
                     hits: List[Chem.Mol],
                     n_cores=1,
                     timeout=600,
                     suffix: str = '',
                     run_plip: bool = True,
                     **settings):
        """
        Redock, but place => replace.
        Ie. score the hits.

        This is not called by ``core_ops``.
        """
        # place themselves for a ∆∆G score
        lab = cls(pdbblock=pdbblock, covalent_resi=None, run_plip=run_plip)
        selfies = pd.DataFrame([dict(name=hit.GetProp('_Name') + '_replaced',
                                     hits=[hit],
                                     smiles=Chem.MolToSmiles(hit)
                                     ) for hit in hits])
        replacements: pd.DataFrame = lab.place(place_input_validator(selfies), n_cores=n_cores, timeout=timeout)
        cls.fix_intxns(replacements)
        replacements['bleached_name'] = replacements['name']
        replacements['name'] = replacements.hit_mols.apply(lambda ms: ms[0].GetProp('_Name'))
        replacements.to_pickle(f'fragmenstein_hit_replacements{suffix}.pkl.gz')
        return replacements
