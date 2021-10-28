PYTHON=/data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python
DATA_ROOT=/data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned
OUTDIR=/data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_Site4


export EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json
export DASK_WORKER_MEMORY="4GB"
$PYTHON  -m examples.protocols_bricsFragmenstein --n_cpus 8 -d $DATA_ROOT -o $OUTDIR -f x0169_0B x290_0B x707_0B

