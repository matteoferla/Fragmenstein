import argparse
from rdkit import Chem
from ..victor import Victor
from ..igor import Igor
from .base import FragmensteinParserBase, _add_common, set_verbose

class FragmensteinParserVictor:

    def _define_victor(self, parser: argparse.ArgumentParser):
        """fragmenstein victor combine -h hit1.mol hit2.mol -t protein.pdb
           fragmenstein victor place -h hit1.mol hit2.mol -s 'CCO'"""
        subparsers = parser.add_subparsers(title='operation', help='combine or place')
        combine_parser = subparsers.add_parser('combine', help='combine')
        _add_common(combine_parser, hits=True, output=True, name=False, template=True)
        combine_parser.set_defaults(func=self.victor_combine)
        place_parser = subparsers.add_parser('place', help='place')
        place_parser.add_argument('-s', '--smiles', required=True)
        _add_common(place_parser, hits=True, output=True, name=True, template=True)
        place_parser.set_defaults(func=self.victor_place)

    def victor_combine(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        Victor.work_path = args.output
        victor = Victor(hits=list(map(Chem.MolFromMolFile, args.hits)), pdb_filename=args.template)
        victor.combine(long_name=args.name)
        return Chem.MolToMolBlock(victor.minimized_mol)

    def victor_place(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        Victor.work_path = args.output
        victor = Victor(hits=list(map(Chem.MolFromMolFile, args.hits)), pdb_filename=args.template)
        victor.place(args.smiles, long_name=args.name)
        return Chem.MolToMolBlock(victor.minimized_mol)
