PYTHON=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python
DATA_ROOT=/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondOthers/fragmenstein/PARP14A_xchemlike/aligned
OUTDIR=/data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/outputPlace_parp14A


export EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json
export PATH=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin:$PATH
export DASK_WORKER_MEMORY="6GB"
export N_CPUS=8
export DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT=200
cmd="$PYTHON  -m examples.protocols_placeFragmenstein -o $OUTDIR -f $OUTDIR/split_input/$1 -i $DATA_ROOT --skip"
#"--filter_out_by_num_inspirational_frags 2"

echo $cmd
exec $cmd
