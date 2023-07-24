"""
This script runs mergers on a provided set of hits,
finds analogues in SmallWorld,
places these and PLIP annotates them.
"""


from typing import List, Dict, Any

import os, re, logging, random, time, argparse
from warnings import warn
from pathlib import Path

import pyrosetta
import pyrosetta_help as ph

import pandas as pd
import pandera.typing as pdt

from smallworld_api import SmallWorld

from rdkit import Chem, rdBase
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

from fragmenstein import Igor, Victor, Laboratory
from fragmenstein.laboratory.validator import place_input_validator

# ------------------------------------------------------

logger = ph.configure_logger()
if logger.handlers:
    logger.handlers[0].setLevel(logging.ERROR)  # logging.WARNING = 30

Igor.init_pyrosetta()

sws = SmallWorld()

chemical_databases: pd.DataFrame = sws.retrieve_databases()

# ------------------------------------------------------

def set_up(output, cutoff, quick, suffix, **kwargs):
    os.makedirs(output, exist_ok=True)
    Victor.work_path = output
    Victor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
    Victor.monster_joining_cutoff = cutoff  # Å
    Victor.quick_reanimation = quick  # for the impatient
    Victor.error_to_catch = Exception  # stop the whole laboratory otherwise
    Victor.enable_stdout(logging.CRITICAL)
    Victor.enable_logfile(os.path.join(output, f'{suffix}.log'), logging.ERROR)


def merge(hits, pdbblock, suffix, n_cores, combination_size, timeout, **settings) -> pd.DataFrame:
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None)
    tick = time.time()
    combinations: pd.DataFrame = lab.combine(hits,
                                             n_cores=n_cores,
                                             timeout=timeout,
                                             combination_size=combination_size)
    with rdBase.BlockLogs():
        combinations['simple_smiles'] = combinations.unminimized_mol.apply(Victor.to_simple_smiles)
    combinations.to_pickle(f'fragmenstein_mergers{suffix}.pkl.gz')
    combinations.to_csv(f'fragmenstein_mergers{suffix}.csv')
    print(tick - time.time())
    return combinations

# ------------------------------------------------------
def get_custom_map(row: pd.Series) -> Dict[str, Dict[int, int]]:
    """
    SmallWorld returns a mapping of the indices.
    This functions maps the indices back to the original hit molecules.

    :param row: SW pd.Series (or dict)
    :return:
    """
    temp = Victor(row.hits, pdb_block=pdbblock)
    temp.monster.positioned_mol = row.unminimized_merger
    temp.minimized_mol = row.minimized_merger
    return temp.migrate_sw_origins(row)

def search(combinations, suffix, sw_dist, sw_length, top_mergers, ranking, sw_db,  **setting) -> pd.DataFrame:
    queries = combinations.sort_values(ranking)\
                          .loc[(combinations.outcome == 'acceptable')]\
                          .drop_duplicates('simple_smiles')\
                          .reset_index()\
                          .head(top_mergers)
    similars = sws.search_many(queries.simple_smiles.to_list(),
                               dist=sw_dist,
                               length=sw_length,
                               db=sw_db,
                               tolerated_exceptions=Exception)
    print(f'Found {len(similars)} analogues')
    # query_index was added clientside to keep track!
    similars['catalogue'] = sw_db
    similars['query_name'] = similars.query_index.map( queries.name.to_dict() )
    similars['hits'] = similars.query_index.map( queries.hit_mols.to_dict() )
    similars['hit_names'] = similars.hits.apply(lambda m: [mm.GetProp('_Name') for mm in m])
    similars['minimized_merger'] = similars.query_index.map( queries.minimized_mol.to_dict() )
    similars['unminimized_merger'] = similars.query_index.map( queries.unminimized_mol.to_dict() )
    similars['name'] = similars['id'] + ':' + similars['query_name']
    similars['smiles'] = similars.hitSmiles.str.split(expand=True)[0]
    similars['custom_map'] = similars.apply(get_custom_map, axis=1)
    similars.to_pickle(f'similars{suffix}.pkl.gz')
    return similars

# ------------------------------------------------------

def place(similars, n_cores, timeout, suffix, **settings) -> pd.DataFrame:
    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=True)
    placements: pd.DataFrame = lab.place(place_input_validator(similars), n_cores=n_cores, timeout=timeout)
    placements.loc[(placements['∆∆G'] > -1) & (placements.outcome == 'acceptable'), 'outcome'] = 'weak'
    placements['unminimized_mol'] = placements.unminimized_mol.fillna(Chem.Mol())
    placements.to_pickle(f'placed{suffix}.pkl.gz')
    placements.to_csv(f'placed{suffix}.csv')
    placements.outcome.value_counts()
    return placements

def write(placements, suffix):
    valids = placements.loc[placements.outcome == 'acceptable'].sort_values('∆∆G')
    valids['combined_name'] = valids['name']
    valids['name'] = valids['enamine_name']
    valids.to_pickle(f'acceptables{suffix}.pkl.gz')
    valids.path.apply(Path.exists).value_counts()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template', help='Template PDB file', required=True)
    parser.add_argument('-i', '--hits', help='Hits SDF file', required=True)
    parser.add_argument('-o', '--output', help='Output folder', default='output')
    parser.add_argument('-r', '--ranking', help='Ranking method', default='∆∆G')
    parser.add_argument('-c', '--cutoff', help='Joining cutoff', default=5)
    parser.add_argument('-q', '--quick', help='Quick reanimation', default=False)
    parser.add_argument('-d', '--sw_dist', help='SmallWorld distance', default=25, type=int)
    parser.add_argument('-l', '--sw_length', help='SmallWorld length', default=50, type=int)
    parser.add_argument('-b', '--sw_databases', help='SmallWorld databases. Accepts multiple e.g. '+\
         '--sw_databases Enamine-SC-Stock-Mar2022.smi.anon Enamine-BB-Stock-Mar2022.smi.anon REAL-Database-22Q1.smi.anon', nargs='+', default=sws.REAL_dataset)
    parser.add_argument('-s', '--suffix', help='Suffix for output files', default='')
    parser.add_argument('-n', '--n_cores', help='Number of cores', default=55, type=int)
    parser.add_argument('-m', '--combination_size', help='Number of hits to combine in one step', default=2, type=int)
    parser.add_argument('-k', '--top_mergers', help='Max number of mergers to followup up on', default=500, type=int)
    parser.add_argument('-e', '--timeout', help='Timeout for each merger', default=240, type=int)
    # load
    settings: Dict[str, Any] = vars(parser.parse_args())
    set_up(**settings)
    with Chem.SDMolSupplier(settings['hits'].strip()) as sd:
        #hitdex: Dict[str, Chem.Mol] = {mol.GetProp('_Name'): mol for mol in sd}
        hits : List[Chem.Mol] = list(sd)
    settings['hits'] = hits
    print(f'N hits: {len(hits)}')
    with open(settings['template'].strip()) as fh:
        pdbblock = fh.read()
    settings['pdbblock'] = pdbblock
    # run
    combinations: pd.DataFrame = merge(**settings)
    uncat_similars: List[pd.DataFrame] = []
    for sw_db in settings['sw_databases']:
        s = search(combinations, sw_db=sw_db, **settings)
        if len(s):
            uncat_similars.append(s)
    assert uncat_similars, 'No analogues were found!'
    similars: pd.DataFrame = pd.concat(uncat_similars, ignore_index=True).drop_duplicates('smiles')
    placements: pd.DataFrame = place(similars, **settings)

