from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    warn('This 3.6+ script **requires** rdkit which ~cannot~ [could not] be pip installed.' +
                              ' To install try either ' +
                              'pip install rdkit-pypi' +
                              'conda install -c conda-forge rdkit or ' +
                              'sudo apt-get/brew install python3-rdkit or visit rdkit documentation.')

if not util.find_spec('pyrosetta'):
    warn('The minimisation part of this code uses pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing. Without it only the classes Monster and Rectifier will work')

if not util.find_spec('pymol2'):
    warn('The module pymol2 is optionally required (conda or apt-get installable).')

setup(
    name='Fragmenstein',
    version='0.6.10',
    packages=find_packages(),
    include_package_data=True,
    package_data={'fragmenstein': ['mpro/data/template.pdb', 'mpro/data/hit_mols/*.mol']},
    install_requires=['pandas', 'numpy', 'rdkit-to-params', 'molecular-rectifier', 'requests'],
    extras_require={'jupyter': ['jupyter']},
    url='https://github.com/matteoferla/Fragmenstein',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Scaffold hopping between bound compounds by stitching them together like a reanimated corpse',
    entry_points={
        'console_scripts': ['fragmenstein=fragmenstein.cli:main'],
    }
)
