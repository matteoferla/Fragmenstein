#!/usr/local/bash
FRAG_ID="x0020"
CONDA_ENVS_DIR=/data/xchem-fragalysis/sanchezg/app/miniconda3/envs/
ACCELLERA_TOOLS=/data/xchem-fragalysis/sanchezg/app/acellera-tools/apps/

ROOT_DATA_DIR=/data/xchem-fragalysis/sanchezg/MD/nsp13/selectedFragments/${FRAG_ID}/
CONDA=/data/xchem-fragalysis/sanchezg/app/miniconda3/condabin/conda
eval "$(${CONDA} shell.bash hook)" &&

cwd=`pwd` &&
cd $ROOT_DATA_DIR &&

echo $CUDA_DEVICE_ORDER
echo $CUDA_VISIBLE_DEVICES

export CUDA_VISIBLE_DEVICES=0

echo "Launching ProteinPreparation" &&
${CONDA_ENVS_DIR}/proteinprepare/bin/python ${ACCELLERA_TOOLS}/ProteinPreparation/app.py -pdb ${ROOT_DATA_DIR}/input/*pdb -outdir ${ROOT_DATA_DIR}/protein_prepared \
                --debug --include-heteroatoms &&

echo "ProteinPreparation done" &&
echo "Launching Parameterize" &&

conda activate paramenv &&
${CONDA_ENVS_DIR}/paramenv/bin/parameterize --no-dihed --charge-type Gasteiger --seed 01 ${ROOT_DATA_DIR}/input/nsp13*mol2 -o ${ROOT_DATA_DIR}/ligand_parametrized &&
#${CONDA_ENVS_DIR}/paramenv/bin/parameterize --no-dihed --charge-type Gasteiger --seed 01 ${ROOT_DATA_DIR}/input/po4.mol2 -o ${ROOT_DATA_DIR}/po4_parametrized&&
conda deactivate &&

echo "Parameterize done" &&
echo "Launching SystemBuilder" &&

conda activate systembuilder &&
${CONDA_ENVS_DIR}/systembuilder/bin/python ${ACCELLERA_TOOLS}/SystemBuilder/app.py -protein ${ROOT_DATA_DIR}/protein_prepared/output.pdb -saltconc 0.1 \
            -ligand ${ROOT_DATA_DIR}/ligand_parametrized/parameters/GAFF2/mol-orig.mol2 \
            -ligfrcmod ${ROOT_DATA_DIR}/ligand_parametrized/parameters/GAFF2/mol.frcmod -outdir ${ROOT_DATA_DIR}/ready_for_MD &&

echo "SystemBuilder done" &&
pwd &&
ls
echo "Launching MDSimulator" &&

conda deactivate &&
conda activate mdsimulator &&

#TODO: add --use-gpu
${CONDA_ENVS_DIR}/mdsimulator/bin/python ${ACCELLERA_TOOLS}/MDSimulator/app.py -inputdir ${ROOT_DATA_DIR}/ready_for_MD/build_1  -constraints protein-ligand  \
            -runtime 50 -equiltime 2 -ligresname UNK &&

cd $cwd
