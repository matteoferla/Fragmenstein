from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Pip and Non pip modules  ----------------------------------------------------------------------------------
requirements = ['pandas', 'numpy', 'rdkit-to-params', 'molecular-rectifier', 'requests']

if not util.find_spec('rdkit'):
    # pypi overwrites the conda version
    requirements.append('rdkit-pypi')

if not util.find_spec('pyrosetta'):
    warn('The minimisation part of this code uses pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing. Without it only the classes Monster and Rectifier will work')

if not util.find_spec('pymol2'):
    warn('The module pymol2 is optionally required (conda or apt-get installable).')

long_description = '''
Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.
<img src="https://github.com/matteoferla/Fragmenstein/blob/master/images/fragmenstein.jpg?raw=true" width="300px">

Documentation in [GitHub](https://github.com/matteoferla/Fragmenstein).

[![colab demo](https://img.shields.io/badge/Run--demo--in--colab-colab_fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab_fragmenstein.ipynb)

![Ox](https://upload.wikimedia.org/wikipedia/en/thumb/2/2f/University_of_Oxford.svg/132px-University_of_Oxford.svg.png)
'''

setup(
    name='Fragmenstein',
    version='0.7.0',
    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
    package_data={'fragmenstein.mpro.data': ['template.pdb'],
                  'fragmenstein.mpro.data.hit_mols': ['*.mol']},
    install_requires=requirements,
    extras_require={'jupyter': ['jupyter']},
    url='https://github.com/matteoferla/Fragmenstein',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    classifiers=[ # https://pypi.org/classifiers/
        'Development Status :: 4 - Beta', # Development Status :: 5 - Production/Stable
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description='Scaffold hopping between bound compounds by stitching them together like a reanimated corpse',
    long_description=long_description,
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': ['fragmenstein=fragmenstein.cli:main'],
    }
)
