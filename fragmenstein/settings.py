default_settings_yaml = """
# These are the default settings for fragmenstein.
# To overide please define $FRAGMENSTEIN_SETTINGS as a yaml file.
# or pass an environment variable _prior_ to import 
# e.g. ff_constraint becomes $FRAGMENSTEIN_FF_CONSTRAINT.
# These will have priority over the defaults.
# Note there are no safeguards against typos.
# For the command line interface, see fragmenstein/_cli_defaults.py

# General settings
work_path: output
monster_average_position: False
monster_throw_on_discard: False
ff_minisation: True

# During the RDKit minisation, how much lee-way to give an atom before it gets penalised.
ff_max_displacement: 0.1

# During the RDKit minisation, how much to penalise an atom that is too far from its ideal position.
ff_constraint: 5.

# During the RDKit minisation, how many iterations to run.
ff_max_iterations: 200

# During the RDKit minisation, use the neighbourhood to constrain the molecule.
ff_use_neighborhood: True
ff_neighborhood: 6.0

# For Wictor, weird things happen if True
ff_minimise_ideal: False

# OpenMM settings
mm_restraint_k: 1000.0
mm_tolerance: 10.0  # mmu.kilocalorie_per_mole / (mmu.nano * mmu.meter)
mm_max_iterations: 0 # 0 is infinite
mm_mobile_radius: 8.0  # mmu.angstrom
"""

import os, yaml
from pathlib import Path
from warnings import warn
from ._cli_defaults import cli_default_settings

if not os.environ.get('FRAGMENSTEIN_SETTINGS', None):
    default_settings = yaml.load(default_settings_yaml, Loader=yaml.FullLoader)
elif Path(os.environ['FRAGMENSTEIN_SETTINGS']).exists():
        with open(os.environ['FRAGMENSTEIN_SETTINGS']) as r:
            default_settings = yaml.load(r, Loader=yaml.FullLoader)
else:
    raise ValueError(f'FRAGMENSTEIN_SETTINGS={os.environ["FRAGMENSTEIN_SETTINGS"]} does not exist as a file')

for variable in os.environ:
    if 'FRAGMENSTEIN_' in variable:
        key = variable.replace('FRAGMENSTEIN_', '').lower()
        if key not in default_settings and key not in cli_default_settings:
            # there should be some kind of check here if a CLI variable is called within Python.
            warn(f'Environment variable {variable} is not a valid setting and will likely be ignored')
        key_type = type(default_settings.get(key, ''))
        default_settings[key] = key_type(os.environ[variable])
