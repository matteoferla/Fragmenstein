#!/usr/bin/env bash

ENVIRON_NAME="Fragmenstein"
conda create -n $ENVIRON_NAME -c rdkit python=3.7 pyrosetta tensorflow-gpu=1.13 rdkit joblib pandas qt requests docopt matplotlib cython nomkl numba tqdm biopython --yes &&

eval "$(conda shell.bash hook)" && conda activate $ENVIRON_NAME &&
pip install dill rdkit_to_params planarity pickle5 molecular_rectifier mrcfile dirsync && conda install -c conda-forge glew glm planarity openbabel lbzip2 --yes &&
conda install -c anaconda netcdf4 --yes && conda install --yes -c conda-forge mesalib &&
git clone https://github.com/schrodinger/pymol-open-source.git && cd pymol-open-source && python setup.py install &&
cd .. && rm -fr pymol-open-source &&
pip install git+https://github.com/pharmai/plip@2d995f872aed9ee829f3cc4f939a21d69a6634b3 --no-deps