import argparse, os, json, itertools, string
import contextlib

from rdkit import Chem
from .._cli_defaults import cli_default_settings
from .base import set_verbose
from typing import List, Any, Dict
from ..laboratory import Laboratory
from ..victor import Victor
import pandas as pd
from rdkit.Chem import PandasTools

class FragmensteinParserPipeline:
    def _define_pipeline(self, parser: argparse.ArgumentParser):
            """fragmenstein pipeline -i hits.sdf -o out.sdf -t template.pdb

            Performs a pipeline run of fragmenstein:
            places the hits against themselves as a reference ("replace"),
            combines the hits,
            searches in SmallWorld for analogues,
            places these,
            and scores the results.
            """
            parser.add_argument('-t', '--template', help='Template PDB file',
                                required=True)
            parser.add_argument('-i', '--input', help='Hits SDF file',
                                required=True)
            parser.add_argument('-o', '--output', help='Output folder',
                                default=cli_default_settings['output'])
            parser.add_argument('-r', '--ranking', help='Ranking method',
                                default=cli_default_settings['ranking'])
            parser.add_argument('-c', '--cutoff', help='Joining cutoff',
                                default=cli_default_settings['cutoff'],
                                type=float)
            parser.add_argument('-q', '--quick',
                                help='Quick reanimation',
                                default=cli_default_settings['quick'],
                                type=bool)
            parser.add_argument('-d', '--sw_dist',
                                help='SmallWorld distance',
                                default=cli_default_settings['sw_dist'],
                                type=int)
            parser.add_argument('-l', '--sw_length', help='SmallWorld length',
                                default=cli_default_settings['sw_length'],
                                type=int)
            parser.add_argument('-b', '--sw_databases', help='SmallWorld databases. Accepts multiple',
                                nargs='+',
                                default=cli_default_settings['sw_databases'])
            parser.add_argument('-s', '--suffix',
                                help='Suffix for output files',
                                default=cli_default_settings['suffix'])
            parser.add_argument('--workfolder',
                                help='Location to put the temp files',
                                default=cli_default_settings['workfolder'])
            parser.add_argument('--victor',
                                help='Which victor to use: Victor, OpenVictor or Wictor',
                                default='Victor')
            parser.add_argument('-n', '--n_cores', help='Number of cores',
                                default=cli_default_settings['n_cores'],
                                type=int)
            parser.add_argument('-m', '--combination_size', help='Number of hits to combine in one step',
                                default=2,
                                type=int)
            parser.add_argument('-k', '--top_mergers', help='Max number of mergers to followup up on',
                                default=cli_default_settings['top_mergers'],
                                type=int)
            parser.add_argument('-e', '--timeout', help='Timeout for each merger',
                                default=cli_default_settings['timeout'],
                                type=int)
            parser.add_argument('-x', '--max_tasks', help='Max number of combinations to try in a batch',
                                default=cli_default_settings['max_tasks'],
                                type=int)
            parser.add_argument('-z', '--blacklist', help='Blacklist file',
                                default=cli_default_settings['blacklist'])
            parser.add_argument('-j', '--weights', help='JSON weights file', default=cli_default_settings['weights'])
            parser.add_argument('-v', '--verbose', action="count", help='verbose')
            parser.set_defaults(func=self.pipeline)
            return parser

    def _pipe_set_up(self, output, cutoff, quick, suffix, **kwargs):
        os.makedirs(output, exist_ok=True)
        Victor.work_path = output
        Victor.monster_throw_on_discard = True  # stop this merger if a fragment cannot be used.
        Victor.monster_joining_cutoff = cutoff  # Å
        Victor.quick_reanimation = quick  # for the impatient
        Victor.error_to_catch = Exception  # stop the whole laboratory otherwise

    def pipeline(self, args: argparse.Namespace):
        """
        Performs a pipeline run of fragmenstein:
        places the hits against themselves as a reference ("replace") to get a baseline score
        combines the hits,
        searches in SmallWorld for analogues,
        places these,
        and scores the results.
        """
        set_verbose(args.verbose)
        settings: Dict[str, Any] = vars(args)
        self._pipe_set_up(**settings)
        with Chem.SDMolSupplier(settings['input'].strip()) as sd:
            # hitdex: Dict[str, Chem.Mol] = {mol.GetProp('_Name'): mol for mol in sd}
            hits: List[Chem.Mol] = list(sd)
        settings['hits'] = hits
        if settings['blacklist']:
            with open(settings['blacklist'].strip()) as fh:
                settings['blacklist'] = [line.strip() for line in fh.readlines()]
        else:
            settings['blacklist'] = cli_default_settings['blacklist']
        print(f'N hits: {len(hits)}')
        with open(settings['template'].strip()) as fh:
            pdbblock = fh.read()
        settings['pdbblock'] = pdbblock
        if settings['weights'] and isinstance(settings['weights'], str):
            with open(settings['weights'].strip()) as fh:
                settings['weights'] = json.load(fh)
        elif settings['weights'] and isinstance(settings['weights'], dict):
            pass  # impossible unless a copy-pasted block.
        else:
            settings['weights'] = cli_default_settings['weights']
        choice = settings.get('victor', 'Victor').lower()
        if choice == 'victor':
            Laboratory.Victor = Victor
        elif choice == 'openvictor':
            from ..openmm.openvictor import OpenVictor
            Laboratory.Victor = OpenVictor
        elif choice == 'wictor':
            from ..faux_victors import Wictor
            Laboratory.Victor = Wictor
        elif choice == 'quicktor':
            from ..faux_victors import Quicktor
            Laboratory.Victor = Quicktor
        else:
            raise ValueError(f'Unknown victor: {choice}')
        Laboratory.Victor.work_path = settings.get('workfolder', cli_default_settings['workfolder'])
        # ## Analyses start here
        # self
        hit_replacements: pd.DataFrame = Laboratory.replace_hits(**settings)
        # run
        max_tasks = settings['max_tasks']
        hitnames = [h.GetProp('_Name') for h in hits]
        all_names = list(map('-'.join, itertools.permutations(hitnames, settings['combination_size'])))
        base_suffix = settings['suffix']
        if max_tasks == 0 or max_tasks > len(all_names):
            all_placements: pd.DataFrame = Laboratory.core_ops(hit_replacements, **settings)
        else:
            all_placements = pd.DataFrame()
            letters = iter(string.ascii_uppercase)
            for i in range(0, len(all_names) + max_tasks, max_tasks):
                settings['suffix'] = base_suffix + next(letters)
                with contextlib.suppress(Exception):
                    placements: pd.DataFrame = Laboratory.core_ops(hit_replacements, **settings)
                    all_placements = pd.concat([all_placements, placements], ignore_index=True)
                settings['blacklist'] += all_names[i:i + max_tasks]
            settings['suffix'] = base_suffix
        Laboratory.correct_weaklings(hit_replacements, all_placements)
        all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
        Laboratory.score(all_placements, hit_replacements, **settings)
        all_placements.to_pickle(f'fragmenstein_placed{base_suffix}.pkl.gz')
        all_placements.to_csv(f'fragmenstein_placed{base_suffix}.csv')
        #PandasTools.WriteSDF(all_placements, f'fragmenstein_placed{base_suffix}.sdf')
        Laboratory.export_sdf(df=all_placements)
