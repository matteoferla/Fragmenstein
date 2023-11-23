from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import os

this_directory = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(this_directory, 'README.md')):
    with open(os.path.join(this_directory, 'README.md'), 'r') as f:
        long_description = f.read()
else:
    long_description = '''
    Scaffold hopping between bound compounds by stitching them together like a reanimated corpse.
    <img src="https://github.com/matteoferla/Fragmenstein/blob/master/images/fragmenstein.jpg?raw=true" width="300px">

    Documentation in [GitHub](https://github.com/matteoferla/Fragmenstein).

    [![colab demo](https://img.shields.io/badge/Run--demo--in--colab-colab_fragmenstein.ipynb-f9ab00?logo=googlecolab)](https://colab.research.google.com/github/matteoferla/Fragmenstein/blob/master/colab_fragmenstein.ipynb)

    ![Ox](https://upload.wikimedia.org/wikipedia/en/thumb/2/2f/University_of_Oxford.svg/132px-University_of_Oxford.svg.png)
    '''

if os.path.exists(os.path.join(this_directory, 'requirements.txt')):
    with open(os.path.join(this_directory, 'requirements.txt'), 'r') as f:
        requirements = [line.split('#')[0].strip() for line in f.readlines()]
        requirements = [line for line in requirements if line]
else:
    requirements = []

# ---------- Pip and Non pip modules  ----------------------------------------------------------------------------------

  # optional

# `pip install xxx` from wheels (the pre-compiled packages), does not run setup.py.
# While sdist does.
# Therefore these warnings will not be shown for pip installations from wheels...


if not util.find_spec('pyrosetta'):
    warn('The minimisation part of this code uses pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing. Without it minimisation via PyRosetta is not possible. '+
         'For openMM usage see OpenVictor.')

if not util.find_spec('pymol2'):
    warn('The module pymol2 is optionally required (conda or apt-get installable).')



setup(
    name='Fragmenstein',
    version='0.13.35',
    description='Merging, linking and placing compounds by stitching them together like a reanimated corpse',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.7',
    packages=find_packages(),
    include_package_data=True,
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
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    entry_points={
        'console_scripts': ['fragmenstein=fragmenstein.cli:main'],
    }
)
