from setuptools import setup
from warnings import warn
from importlib import util

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    raise ModuleNotFoundError('This 3.6+ script **requires** rdkit which cannot be pip installed.' +
                              ' To install try either ' +
                              'conda install -c conda-forge rdkit or ' +
                              'sudo apt-get/brew install python3-rdkit or visit rdkit documentation.')

if not util.find_spec('pyrosetta'):
    warn('The minimisation part of this code uses pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing. Without it only the classes Monster and Rectifier will work')

if not util.find_spec('pymol2'):
    warn('The module pymol2 is optionally required (conda or apt-get installable).')

setup(
    name='Fragmenstein',
    version='0.5',
    packages=['fragmenstein'],
    install_requires=['numpy'],
    extras_require={'minimization': ['rdkit_to_params'],
                    'jupyter': ['jupyter']},
    url='https://github.com/matteoferla/Fragmenstein',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Scaffold hopping between bound compounds by stitching them together like a reanimated corpse',
    entry_points={
        'console_scripts': ['fragmenstein=fragmenstein.cli:main'],
    }
)
