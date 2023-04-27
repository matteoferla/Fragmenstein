import argparse, logging
from ..victor import Victor
from typing import Optional, List


class FragmensteinParserBase:
    """
    See main module docstring for usage.
    """

    subcommand_help = '''Actions: monster, victor, laboratory, utils
    '''
    monster_help = '''monster'''
    laboratory_help = '''laboratory'''
    victor_help = '''victor'''
    utils_help = '''utils'''

    def __init__(self):
        self.parser = argparse.ArgumentParser(description=self.__doc__)
        subparsers = self.parser.add_subparsers(title='subcommands', help=self.subcommand_help)
        monster_parser = subparsers.add_parser('monster', help=self.monster_help)
        self._define_monster(monster_parser)
        laboratory_parser = subparsers.add_parser('laboratory', help=self.laboratory_help)
        self._define_laboratory(laboratory_parser)
        victor_parser = subparsers.add_parser('victor', help=self.victor_help)
        self._define_victor(victor_parser)
        utils_parser = subparsers.add_parser('utils', help=self.utils_help)
        self._define_utils(utils_parser)

    def __call__(self, cli_override: Optional[List[str]] = None):
        if cli_override:
            args = self.parser.parse_args(cli_override)
        else:
            args = self.parser.parse_args()
        args.func(args)

    def _define_laboratory(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_victor(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_monster(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')

    def _define_utils(self, parser: argparse.ArgumentParser):
        raise NotImplementedError('virtual method')


def set_verbose(count):
    log_map = {None: logging.FATAL,
               0: logging.CRITICAL,
               1: logging.ERROR,
               2: logging.WARNING,
               3: logging.INFO,
               4: logging.DEBUG}
    Victor.enable_stdout(log_map.get(count, logging.DEBUG))
    Victor.capture_logs()

def _add_common(parser: argparse.ArgumentParser, **settings):
    if settings.get('verbose', True):
        parser.add_argument('-v', '--verbose', action="count", help='verbose')
    if settings.get('hits', False):
        parser.add_argument('-i', '--hits', nargs='+', required=True, help='hit mol files')
    if settings.get('output', False):
        parser.add_argument('-o', '--output', default='.', help='output root folder')
    if settings.get('name', False):
        parser.add_argument('-n', '--name', default='fragmenstein', help='output name of molecule')
    if settings.get('template', False):
        parser.add_argument('-t', '--template', required=True, help='template PDB file')
