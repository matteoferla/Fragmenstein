import argparse
import pandas as pd
from .base import FragmensteinParserBase, _add_common, set_verbose
from rdkit import Chem
from rdkit.Chem import PandasTools
from ..laboratory import Laboratory
from ..igor import Igor
from typing import List, Tuple

class FragmensteinParserLaboratory:

    def common(self, parser: argparse.ArgumentParser):
        # hits in _add_common is for individual mol files
        _add_common(parser, hits=False, output=True, template=True)
        parser.add_argument('-i', '--input', required=True, help='input sdf file')
        parser.add_argument('-d', '--out-table', default='output.csv', help='table output file')
        parser.add_argument('-s', '--sdf-outfile', default='output.sdf', help='sdf output file')
        parser.add_argument('-c', '--cores', default=1, type=int, help='number of cores to use')
        parser.add_argument('-p', '--run-plip', default=False, type=bool, help='Run PLIP?')
        parser.add_argument('--victor',
                            help='Which victor to use: Victor, OpenVictor or Wictor',
                            default='Victor')

    def _define_laboratory(self, parser: argparse.ArgumentParser):
        """fragmenstein laboratory combine -i hits.sdf -o out.sdf
           fragmenstein laboratory place -i hits.sdf -t table.csv -o out.sdf"""
        subparsers = parser.add_subparsers(title='operation', help='combine or place')
        # combine
        combine_parser = subparsers.add_parser('combine', help='combine')
        self.common(combine_parser)
        combine_parser.set_defaults(func=self.lab_combine)
        # place
        place_parser = subparsers.add_parser('place', help='place')
        self.common(place_parser)
        place_parser.add_argument('-f', '--in-table', required=True,
                                  help='CSV table input file, ' +\
                                       'requires `name`, `smiles` and space-separated-`hit_names`')
        place_parser.set_defaults(func=self.lab_place)

    def gather(self, args: argparse.Namespace) -> Tuple[Laboratory, List[Chem.Mol]]:
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        # victor or wictor? copypasted from pipeline....
        choice = args.victor.lower()
        if choice == 'victor':
            from .victor import Victor
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
        # laboratory
        lab = Laboratory(pdbblock=open(args.template).read(), run_plip=bool(args.run_plip))
        with Chem.SDMolSupplier(args.input) as suppl:
            mols = list(suppl)
        return lab, mols

    def lab_combine(self, args: argparse.Namespace) -> str:
        lab, mols = self.gather(args)
        combos: pd.DataFrame = lab.combine(mols, n_cores=args.cores)
        self.write(args, combos)

    def write(self, args: argparse.Namespace, df: pd.DataFrame):
        props = ['smiles', 'error', '∆∆G', '∆G_bound', '∆G_unbound', 'comRMSD',
                 'N_constrained_atoms', 'N_unconstrained_atoms', 'runtime', 'regarded',
                 'disregarded', 'LE', 'outcome',
                 'percent_hybrid']
        df[props].to_csv(args.out_table)
        PandasTools.WriteSDF(df, args.sdf_outfile,
                             'minimized_mol',
                             properties=props)

    def fix_columns(self, df: pd.DataFrame, preferred_name: str, options: List[str]):
        if preferred_name not in df.columns:
            for name_col in options:
                for variant_fx in [str.lower, str.upper, str.title]:
                    alt = variant_fx(name_col)  # noqa
                    if alt in df.columns:
                        df.rename(columns={alt: preferred_name}, inplace=True)
                        break
            else:
                raise ValueError(f'No name column `{preferred_name}` (or the fallbacks {options}) in {df.columns}')

    def lab_place(self, args: argparse.Namespace) -> str:
        lab, hits = self.gather(args)
        hitdex = {hit.GetProp('_Name'): hit for hit in hits}
        # fix the input table
        intable = pd.read_csv(args.in_table)
        # name
        self.fix_columns(intable, 'name', ['name', 'long name'])
        self.fix_columns(intable, 'smiles', ['smiles', 'smile'])
        self.fix_columns(intable, 'hit_names', ['hit_names', 'hit names', 'hits'])
        try:
            intable['hits'] = intable['hit_names'].apply(lambda x: [hitdex[n] for n in x.split()])
        except KeyError as error:
            raise KeyError(f'Could not find {error} in the input sdf file ({hitdex.keys()})')
        # run
        placements = lab.place(intable, n_cores=args.cores)
        self.write(args, placements)



