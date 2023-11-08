from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd
from rdkit.Chem import Descriptors

from typing import List, Dict, Any, Optional
import operator, os, re, logging, random, time, argparse, string, itertools, json, contextlib, requests
from warnings import warn
from ..monster import Monster
import pandas as pd
import pandera.typing as pdt
from pandarallel import pandarallel
from smallworld_api import SmallWorld


from rdkit import Chem, rdBase, DataStructs
from rdkit.Chem import AllChem, Draw, PandasTools, BRICS
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFingerprintGenerator as rdfpg
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams, FilterCatalog, FilterCatalogEntry
from .._cli_defaults import cli_default_settings

# ----- Scoring -----------------------------------
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog(params)

def get_pains(mol) -> List[str]:
    with contextlib.suppress(Exception):
        entry: FilterCatalogEntry
        if not isinstance(mol, Chem.Mol) or mol.GetNumHeavyAtoms() == 0:
            return []
        AllChem.SanitizeMol(mol)
        return [entry.GetDescription() for entry in catalog.GetMatches(mol)]

class GetRowSimilarity:
    def __init__(self, hits):
        self.fpgen = rdfpg.GetRDKitFPGenerator()
        self.hit2fp = {h.GetProp('_Name'): self.fpgen.GetFingerprint(h) for h in hits}

    def __call__(self, row: pd.Series):
        with contextlib.suppress(cli_default_settings['supressed_exceptions']):
            if not isinstance(row.minimized_mol, Chem.Mol):
                return float('nan')
            elif isinstance(row.hit_names, str):
                hit_names = row.hit_names.split(',')
            elif isinstance(row.hit_names, list):
                hit_names = row.hit_names
            else:
                return float('nan')
            fp = self.fpgen.GetFingerprint(AllChem.RemoveHs(row.minimized_mol))
            return max([DataStructs.TanimotoSimilarity(fp, self.hit2fp[name]) for name in hit_names])
        return float('nan')


class HitIntxnTallier:

    def __init__(self, hit_replacements):
        self.slim_hits = self.slim_down(hit_replacements)

    def slim_down(self, hit_replacements):
        # the bleaching was fixed cf bleached_name
        # hit_replacements['new_name'] = hit_replacements.name.str.replace('-', '_')
        # undoing bleaching
        hit_replacements['new_name'] = hit_replacements.hit_mols.apply(lambda ms: ms[0].GetProp('_Name'))
        columns = [c for c in hit_replacements.columns if isinstance(c, tuple)]
        return hit_replacements.set_index('new_name')[columns].fillna(0).copy()

    def __call__(self, row: pd.Series):
        with contextlib.suppress(cli_default_settings['supressed_exceptions']):
            if not isinstance(row.minimized_mol, Chem.Mol) or isinstance(row.hit_names, float):
                return float('nan'), float('nan')
            present_tally = 0
            absent_tally = 0
            for hit_name in list(row.hit_names):
                if hit_name not in self.slim_hits.index:
                    raise Exception('Name' + hit_name)
                hit_row = self.slim_hits.loc[hit_name]
                for intxn_name, hit_value in hit_row.items():
                    if not hit_value:
                        continue
                    elif intxn_name not in row.index:
                        absent_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
                    elif row[intxn_name]:
                        absent_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
                    else:
                        present_tally += 1 if intxn_name[0] != 'hydroph_interaction' else 0.5
            return present_tally, absent_tally
        return float('nan'), float('nan')


class UniquenessMeter:
    def __init__(self, tallies, intxn_names, k=0.5):
        self.tallies = tallies
        self.intxn_names = intxn_names
        self.k = k

    def __call__(self, row):
        with contextlib.suppress(cli_default_settings['supressed_exceptions']):
            return sum([(row[name] / self.tallies[name]) ** self.k for name in self.intxn_names if
                        row[name] and self.tallies[name]])
        return float('nan')

    def tally_interactions(self, row):
        return sum([row[c] if self.intxn_names[0] != 'hydroph_interaction' else row[c] * 0.5 for c in self.intxn_names])


class PenaltyMeter:
    def __init__(self, weights, nan_penalty=10):
        self.weights = weights
        self.nan_penalty = nan_penalty

    def __call__(self, row):
        with contextlib.suppress(cli_default_settings['supressed_exceptions']):
            penalty = 0
            if row.outcome != 'acceptable':
                return float('inf')
            for col, w in self.weights.items():
                if col not in row.index:
                    warn(f'{col} column is missing from df')
                    continue
                penalty += row[col] * w if str(row[col]) != 'nan' else self.nan_penalty
            return penalty
        return float('nan')


def butina_cluster(mol_list, cutoff=0.35):
    # https://github.com/PatWalters/workshop/blob/master/clustering/taylor_butina.ipynb
    fp_list = [rdmd.GetMorganFingerprintAsBitVect(AllChem.RemoveAllHs(m), 3, nBits=2048) for m in mol_list]
    dists = []
    nfps = len(fp_list)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dists.extend([1 - x for x in sims])
    mol_clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    cluster_id_list = [0] * nfps
    for idx, cluster in enumerate(mol_clusters, 1):
        for member in cluster:
            cluster_id_list[member] = idx
    return cluster_id_list

def UFF_Gibbs(mol):
    # free energy cost of bound conformer
    if not isinstance(mol, Chem.Mol) or mol.GetNumHeavyAtoms() == 0:
        return float('nan')
    with contextlib.suppress(cli_default_settings['supressed_exceptions']):
        AllChem.SanitizeMol(mol)
        # this is actually UFF
        copy = Chem.Mol(mol)
        return Monster.MMFF_score(None, mol, True)
    return float('nan')

class LabScore:

    @classmethod
    def score(cls,
              placements: pd.DataFrame,
              hit_replacements: pd.DataFrame,
              weights: dict,
              **settings):
        """
        This is very much a method for the CLI.
        A real Pythonic usage would be to address the individual components.
        """
        # tanimoto
        hits: List[Chem.Mol] = hit_replacements.hit_mols.apply(operator.itemgetter(0)).to_list()
        get_similarity = GetRowSimilarity(hits)
        placements['max_hit_Tanimoto'] = placements.apply(get_similarity, axis=1)
        # properties
        m = placements.minimized_mol.apply(lambda m: m if isinstance(m, Chem.Mol) else Chem.Mol())
        # macrocyclics... yuck.
        placements['largest_ring'] = m.apply(lambda mol: max([0] + list(map(len, mol.GetRingInfo().AtomRings()))))
        # interactions
        cls.fix_intxns(placements)
        tally_hit_intxns = HitIntxnTallier(hit_replacements)
        hit_checks = placements.apply(tally_hit_intxns, axis=1)
        placements['N_interactions_kept'] = hit_checks.apply(operator.itemgetter(0))  # .fillna(0).astype(int)
        placements['N_interactions_lost'] = hit_checks.apply(operator.itemgetter(1))  # .fillna(99).astype(int)
        intxn_names = [c for c in placements.columns if isinstance(c, tuple)]
        tallies = placements[intxn_names].sum()
        ratioed = UniquenessMeter(tallies, intxn_names, k=0.5)
        placements['interaction_uniqueness_metric'] = placements.apply(ratioed, axis=1)
        placements['N_interactions'] = placements.apply(ratioed.tally_interactions, axis=1)
        placements['PAINSes'] = placements.minimized_mol.apply(get_pains)
        placements['N_PAINS'] = placements.PAINSes.apply(len)
        placements['UFF_Gibbs'] = placements.minimized_mol.apply(UFF_Gibbs)
        placements['strain_per_HA'] = placements.UFF_Gibbs / (placements.N_HA + 0.0001)
        penalize = PenaltyMeter(weights)
        placements['ad_hoc_penalty'] = placements.apply(penalize, axis=1)
        with contextlib.suppress(cli_default_settings['supressed_exceptions']):
            placements['cluster'] = butina_cluster(m.to_list())

    @staticmethod
    def export_sdf(df: pd.DataFrame,
                   penalty_col:str = 'ad_hoc_penalty',
                   filename: str = 'fragmenstein.sdf',
                   target_name: str = 'DOESNT_MATTER_UNLESS_YOU_WANT_TO_REMOVE_IT'):
        def fix(mol: Chem.Mol) -> None:
            assert isinstance(mol, Chem.Mol)
            assert mol.GetNumAtoms()
            mol.ClearComputedProps()
            for name in mol.GetPropNames():
                mol.ClearProp(name)

        with pd.option_context('mode.chained_assignment', None):
            df = df.loc[df.outcome == 'acceptable'] \
                .sort_values(penalty_col) \
                .rename(columns={c: ':'.join(map(str, c)) for c in df.columns if isinstance(c, tuple)}) \
                .reset_index() \
                .copy()
            # list of str to str w/ comma-separator
            df['ref_mols'] = df.hit_names.apply(lambda l: ','.join([v.replace(f'{target_name}-', '') for v in l]))
            df['washed_mol'] = df.minimized_mol.apply(fix)
            df['name'] = df['name'].apply(lambda v: v.split('-D68EV3CPROA')[0])
            # non str/float/ints
            not_okay = ('name', 'minimized_mol', 'ref_mols', 'washed_mol', 'mode',
                        'runtime', 'error', 'outcome',
                        'smiles',
                        'regarded', 'disregarded', 'hit_names',
                        '∆G_bound',
                        '∆G_unbound',
                        'unmin_binary',
                        'min_binary',
                        'hit_binaries',
                        'minimized_mol',
                        'hit_mols', 'unminimized_mol', 'hit_names')
        good_columns = df.columns[~df.map(lambda x: not isinstance(x, (float, str))).any()]
        extras: List[str] = [c for c in df.columns if c in good_columns and not c in not_okay]
        PandasTools.WriteSDF(df, out=filename, properties=extras, molColName='minimized_mol')