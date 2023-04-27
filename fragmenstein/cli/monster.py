from .base import FragmensteinParserBase, _add_common, set_verbose
import argparse
from fragmenstein import Monster
from rdkit import Chem

class FragmensteinParserMonster:

    def _define_monster(self, parser: argparse.ArgumentParser):
        """fragmenstein monster combine -i hit1.mol hit2.mol
           fragmenstein monster place -i hit1.mol hit2.mol -s 'CCO'"""
        subparsers = parser.add_subparsers(title='operation', help='combine or place')
        combine_parser = subparsers.add_parser('combine', help='combine')
        _add_common(combine_parser, hits=True, output=False, name=False)
        combine_parser.set_defaults(func=self.monster_combine)
        place_parser = subparsers.add_parser('place', help='place')
        place_parser.add_argument('-s', '--smiles', required=True)
        _add_common(place_parser, hits=True, output=False, name=True)
        place_parser.set_defaults(func=self.monster_place)

    def monster_combine(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        monsta = Monster(hits=list(map(Chem.MolFromMolFile, args.hits)))
        monsta.combine()
        return Chem.MolToMolBlock(monsta.positioned_mol)

    def monster_place(self, args: argparse.Namespace) -> str:
        set_verbose(args.verbose)
        monsta = Monster(hits=list(map(Chem.MolFromMolFile, args.hits)))
        monsta.place_smiles(args.smiles, long_name=args.name)
        return Chem.MolToMolBlock(monsta.positioned_mol)
