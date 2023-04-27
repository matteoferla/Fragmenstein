#!/usr/bin/env python3

from ..version import __version__
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

from .monster import FragmensteinParserMonster
from .parser import FragmensteinParser
from .base import FragmensteinParserBase

if __name__ == '__main__':
    FragmensteinParser.__doc__ = __doc__
    FragmensteinParser()
