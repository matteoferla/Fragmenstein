import argparse
from .base import set_verbose
from pathlib import Path
from ..igor import Igor, pyrosetta
from collections import Counter

class FragmensteinParserUtils:

    def _define_utils(self, parser: argparse.ArgumentParser):
        subparsers = parser.add_subparsers(title='operation', help='minimize or ...')
        minimize_parser = subparsers.add_parser('minimize',
                                               help='prepare (minimise) the structure and get rid of waters etc')
        minimize_parser.set_defaults(func=self.igor_minimize)
        minimize_parser.add_argument('-t', '--template', required=True, help='Template PDB file')
        minimize_parser.add_argument('-o', '--output', default='', help='output PDB folder')
        minimize_parser.add_argument('-v', '--verbose', action="count", help='verbose')
        minimize_parser.add_argument('-ed', '--electron_density', default='', help='electron density map')
        minimize_parser.add_argument('-cw', '--constraint_weight', default=5, help='constraint weight')
        minimize_parser.add_argument('-c', '--cycles', default=15, help='number of cycles')
        minimize_parser.add_argument('-cf', '--constraint_file', default='', help='constraint file')
        minimize_parser.add_argument('-first', '--first_chain_only', default=False, help='only keep first chain')

    def igor_minimize(self, args: argparse.Namespace):
        set_verbose(args.verbose)
        Igor.init_pyrosetta()
        pyrosetta.rosetta.basic.options.set_boolean_option('run:ignore_zero_occupancy', False)
        pyrosetta.rosetta.basic.options.set_boolean_option('in:auto_setup_metals', True)
        # ## Prep
        assert Path(args.template).exists(), f'{args.template} does not exist'
        if args.output == '':
            args.output = args.template.replace('.pdb', '_min.pdb')
        assert not Path(args.outfile).exists(), f'{args.outfile} already exists'
        # ## Read pose
        pose: pyrosetta.Pose = pyrosetta.pose_from_file(args.template)
        con_tally = Counter(
            [type(con).__name__ for con in pose.constraint_set().get_all_constraints()]).most_common()
        print('chains', pyrosetta.rosetta.core.conformation.pose.conf2pdb_chain(pose))
        print('sequence', pose.sequence())
        print('Constraints present', con_tally)
        if args.first_chain_only:
            original: pyrosetta.Pose = pose.split_by_chain(1)
        if args.constraint_file:
            Igor.add_constraint_file(pose,
                                     constraint_file=args.constraint_file)
        # ## Relax
        if args.electron_density == '':
            relaxed: pyrosetta.Pose = Igor.relax_with_ED(pose,
                                                           ccp4_file=args.ccp4_file,
                                                           constraint_weights=(args.constraint_weight,),
                                                           cycles=args.cycles)
        else:
            relaxed: pyrosetta.Pose = Igor.tethered_relax(pose,
                                                          constraint_weight=args.constraint_weight,
                                                          cycles=args.cycles)
        # ## Dump
        relaxed.dump_pdb(args.outfile)
