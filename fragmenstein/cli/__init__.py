#!/usr/bin/env python3

__version__ = '0.1'
__doc__ = """
Command line interface to Fragmenstein.
The the first argument to the command is one of the three option:

## extract
Given a target folder (`-i <folder_name>`) 
extract the mol files (with covalent atom if present) in an output folder (`-o <folder_name>`) 

## merge


## position

NB. Please no spaces in the filenames. If spaces are present quote/double-quote the fle path.
""".strip()

from typing import *
import argparse


# from fragmenstein import Victor

class FragmensteinParser:
    ___doc__ = __doc__  # this works (it is in a function it doesn't as the global cmd is called for)

    # ===================== Actions ====================================================================================

    def extract(self, input_folder: str, output_folder: str):
        pass

    def merge(self, target_filenames: List):
        pass

    def place(self, target_filenames: List, smiles: str):
        pass

    def set_verbosity(verbose: int):
        # set logging stdout
        # logging.ERROR, logging.CRITICAL --> 0 ?
        # logging.WARNING --> 1 ?
        # logging.INFO --> 2
        # logging.DEBUG --> 3
        pass

    # ===================== Arguments ==================================================================================

    def __init__(self, cli_override: Optional[List[str]] = None):
        self.parser = self.define_parser()
        # ---- parse ---------------------------------------------------------------------------------------------------
        if cli_override:  # testing basically
            self.args = self.parser.parse_args(cli_override)
        else:  # normal way
            self.args = self.parser.parse_args()
        # ---- verify --------------------------------------------------------------------------------------------------
        self.verify()
        # ---- route ---------------------------------------------------------------------------------------------------
        self.route()

    @classmethod
    def mock(cls):
        self = cls.__new__(cls)
        self.parser = self.define_parser()
        self.args = property(lambda: exec("raise NotImplementedError"))

    def define_parser(self):
        parser = argparse.ArgumentParser(description=self.__doc__)
        parser.add_argument('mode', type=str, help='Three modes are accepted: extract, merge, position',
                            choices=['extract', 'merge', 'position'])
        parser.add_argument('-v', '--verbose', action="count", help='verbose')
        # ------- mode: extract ----------------------------------------------------------------------------------------
        parser.add_argument('-i', '--input',
                            type=str, nargs='+', help='mode: extract. input folder')
        parser.add_argument('-o', '--output',
                            type=str, nargs='+',
                            help='mode: extract. output folder')
        # ------- mode: merge and place --------------------------------------------------------------------------------
        parser.add_argument('-t', '--targets',
                            type=str, nargs='+',
                            help='mode: merge and place. input target filenames')
        parser.add_argument('-n', '--name',
                            type=str, nargs=1,
                            help='mode: merge and place. output compound name')
        # ------- mode: place ------------------------------------------------------------------------------------------
        parser.add_argument('-s', '--smiles', type=str, nargs='+', help='mode: place. SMILE-String to be placed')
        return parser

    def verify(self):
        if self.args.mode == 'extract' and not (self.args.input and self.args.output):
            self.parser.error('mode extract requires --input <folder> and --output <folder>')
        elif self.args.mode == 'extract' and not self.args.targets:
            self.parser.error('mode extract requires --targets <file names>')
        elif self.args.mode == 'place' and not (self.args.targets and self.args.smiles):
            self.parser.error('mode place requires --targets <file names> and --smiles <string>')
        else:
            pass  # no issues.

    def route(self):
        self.set_verbosity(self.args.verbose)
        if self.args.mode == 'extract':
            self.extract(self.args.input, self.args.output)
        elif self.args.mode == 'merge':
            self.merge(self.args.targets)
        elif self.args.mode == 'place':
            self.place(self.args.targets, self.args.smiles)
        else:
            raise SyntaxError('Mode choice impossible. Check declaration')


# ======================================================================================================================

def main():
    FragmensteinParser()


if __name__ == '__main__':
    main()
