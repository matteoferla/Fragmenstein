import argparse
import pandas as pd
from .base import FragmensteinParserBase, _add_common, set_verbose
from rdkit import Chem
from rdkit.Chem import PandasTools
from ..laboratory import Laboratory
from ..igor import Igor

class FragmensteinParserLaboratory:

    def _define_laboratory(self, parser: argparse.ArgumentParser):
        """fragmenstein laboratory combine -i hits.sdf -o out.sdf
           fragmenstein laboratory place -i hits.sdf -t table.csv -o out.sdf"""
        subparsers = parser.add_subparsers(title='operation', help='combine or place')
        combine_parser = subparsers.add_parser('combine', help='combine')
        _add_common(combine_parser, hits=False, output=True, template=True)
        parser.add_argument('-i', '--input', required=True, help='input sdf file')
        parser.add_argument('-d', '--table', default='output.csv', help='table output file')
        parser.add_argument('-s', '--sdf_outfile', default='output.sdf', help='sdf output file')
        parser.add_argument('-c', '--cores', default=1, type=int, help='number of cores to use')
        combine_parser.set_defaults(func=self.lab_combine)
        place_parser = subparsers.add_parser('place', help='place')
        place_parser.set_defaults(func=self.lab_place)

    def lab_combine(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        lab = Laboratory(pdbblock=open(args.template).read(), )
        with Chem.SDMolSupplier(args.input) as suppl:
            mols = list(suppl)
        combos: pd.DataFrame = lab.combine(mols, n_cores=args.cores)
        props = ['smiles', 'error','∆∆G', '∆G_bound', '∆G_unbound', 'comRMSD',
       'N_constrained_atoms', 'N_unconstrained_atoms', 'runtime', 'regarded',
       'disregarded', 'LE', 'outcome',
       'percent_hybrid']
        combos[props].to_csv(args.table)
        PandasTools.WriteSDF(combos, args.sdf_outfile,
                             'minimized_mol',
                             properties=props)

    def lab_place(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        lab = Laboratory(pdbblock=open(args.template).read(), )
        raise NotImplementedError('MF has not thought hard enough how this could be done')
