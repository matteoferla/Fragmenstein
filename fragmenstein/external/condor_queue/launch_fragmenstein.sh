PYTHON=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python
DATA_ROOT=/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned
OUTDIR=/data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output


export EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json
export DASK_WORKER_MEMORY="4GB"
$PYTHON  -m examples.protocols_bricsFragmenstein --n_cpus 8 -d $DATA_ROOT -o $OUTDIR -f x0176_0B x0246_0B x0438_0B



