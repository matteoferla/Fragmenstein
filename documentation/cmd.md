# Command line interface

The strength of Fragmenstein is as a python module, but there is a command line interface.
The terminal command `fragmenstein` is the main entry point, and has the following subcommands:

* `utils`: A few utilities, such as minimising a PDB file.
* `monster`: stitching (by combining/placing) compounds regardless of protein, yielding monster molecules.
* `victor`: stitching (by combining/placing) compounds in a protein
* `laboratory`: combinatorial operations
* `pipeline`: The full pipeline, which does the whole thing.

In addition to the command line arguments there are environment variables that can be set, see below.

These in turn have further subcommands:

## Utils

Example:

```bash
fragmenstein utils minimize -t apo_protein.pdb -o minimised.pdb;
fragmenstein utils minimize -t apo_protein.pdb -o minimised.pdb -first;  # only keep first chain
fragmenstein utils minimize -t apo_protein.pdb -o minimised.pdb 0ed map.ccp4 -cw 10 -c 15;  # electron density map
```

Command:

    usage: fragmenstein utils minimize [-h] -t TEMPLATE [-o OUTPUT] [-v] [-ed ELECTRON_DENSITY] [-cw CONSTRAINT_WEIGHT] [-c CYCLES] [-cf CONSTRAINT_FILE] [-first FIRST_CHAIN_ONLY]
    
    options:
      -h, --help            show this help message and exit
      -t TEMPLATE, --template TEMPLATE
                            Template PDB file
      -o OUTPUT, --output OUTPUT
                            output PDB folder
      -v, --verbose         verbose
      -ed ELECTRON_DENSITY, --electron_density ELECTRON_DENSITY
                            electron density map
      -cw CONSTRAINT_WEIGHT, --constraint_weight CONSTRAINT_WEIGHT
                            constraint weight
      -c CYCLES, --cycles CYCLES
                            number of cycles
      -cf CONSTRAINT_FILE, --constraint_file CONSTRAINT_FILE
                            constraint file
      -first FIRST_CHAIN_ONLY, --first_chain_only FIRST_CHAIN_ONLY
                            only keep first chain

## Monster combine

See [workings.md](workings.md) for details.

Example:

```bash
fragmenstein monster combine -i mol1.mol mol2.mol >> combo.mol;
```

Usage:

    usage: fragmenstein monster combine [-h] [-v] -i HITS [HITS ...]
    
    options:
      -h, --help            show this help message and exit
      -v, --verbose         verbose
      -i HITS [HITS ...], --hits HITS [HITS ...]
                            hit mol files

## Monster place

See [workings.md](workings.md) for details.

Example:

```bash
fragmenstein monster place -i mol1.mol mol2.mol >> combo.mol;
```

Usage:

    usage: fragmenstein monster place [-h] -s SMILES [-v] -i HITS [HITS ...] [-n NAME]
    
    options:
      -h, --help            show this help message and exit
      -s SMILES, --smiles SMILES
      -v, --verbose         verbose
      -i HITS [HITS ...], --hits HITS [HITS ...]
                            hit mol files
      -n NAME, --name NAME  output name of molecule

## Victor combine

    
See [workings.md](workings.md) and [victor.md](victor.md) for details.

Example:

```bash
fragmenstein victor combine -i mol1.mol mol2.mol -t protein.pdb -o output_folder >> combo.mol;
```

Usage:

    usage: fragmenstein victor combine [-h] [-v] -i HITS [HITS ...] [-o OUTPUT] -t TEMPLATE
    
    options:
      -h, --help            show this help message and exit
      -v, --verbose         verbose
      -i HITS [HITS ...], --hits HITS [HITS ...]
                            hit mol files
      -o OUTPUT, --output OUTPUT
                            output root folder
      -t TEMPLATE, --template TEMPLATE
                            template PDB file

## Victor place

    
See [workings.md](workings.md) and [victor.md](victor.md) for details.

Example:

```bash
fragmenstein victor place -i mol1.mol mol2.mol -s 'Cn1cnc2c1c(=O)n(C)c(=O)n2C' -t protein.pdb -o output_folder >> placed.mol;
```

Usage:

    usage: fragmenstein victor place [-h] -s SMILES [-v] -i HITS [HITS ...] [-o OUTPUT] [-n NAME] -t TEMPLATE
    
    options:
      -h, --help            show this help message and exit
      -s SMILES, --smiles SMILES
      -v, --verbose         verbose
      -i HITS [HITS ...], --hits HITS [HITS ...]
                            hit mol files
      -o OUTPUT, --output OUTPUT
                            output root folder
      -n NAME, --name NAME  output name of molecule
      -t TEMPLATE, --template TEMPLATE
                            template PDB file


## Laboratory and pipeline subcommands

See [workings.md](workings.md) and [pipeline.md](pipeline.md) for details.

For laboratory and pipeline subcommands, the number of cores can be set.
When not specified, the number of core of the machine is used, this is NOT the number of available cores.
If more cores are set than available, the timeout will affect the run,
so when running on a shared node in slurm, then set `-c $SLURM_CPUS_ON_NODE`.

### Laboratory combine
Example: 

```bash
fragmenstein laboratory combine -i hits.sdf ....;
```

Usage:

    usage: fragmenstein laboratory combine [-h] [-v] [-o OUTPUT] -t TEMPLATE -i INPUT [-d OUT_TABLE] [-s SDF_OUTFILE] [-c CORES] [-p RUN_PLIP] [--victor VICTOR]
    
    options:
      -h, --help            show this help message and exit
      -v, --verbose         verbose
      -o OUTPUT, --output OUTPUT
                            output root folder
      -t TEMPLATE, --template TEMPLATE
                            template PDB file
      -i INPUT, --input INPUT
                            input sdf file
      -d OUT_TABLE, --out-table OUT_TABLE
                            table output file
      -s SDF_OUTFILE, --sdf-outfile SDF_OUTFILE
                            sdf output file
      -c CORES, --cores CORES
                            number of cores to use
      -p RUN_PLIP, --run-plip RUN_PLIP
                            Run PLIP?
      --victor VICTOR       Which victor to use: Victor, OpenVictor or Wictor

### Laboratory place

Example:

```bash
...
```

Usage:

    usage: fragmenstein laboratory place [-h] [-v] [-o OUTPUT] -t TEMPLATE -i INPUT [-d OUT_TABLE] [-s SDF_OUTFILE] [-c CORES] [-p RUN_PLIP] [--victor VICTOR] -f IN_TABLE
    
    options:
      -h, --help            show this help message and exit
      -v, --verbose         verbose
      -o OUTPUT, --output OUTPUT
                            output root folder
      -t TEMPLATE, --template TEMPLATE
                            template PDB file
      -i INPUT, --input INPUT
                            input sdf file
      -d OUT_TABLE, --out-table OUT_TABLE
                            table output file
      -s SDF_OUTFILE, --sdf-outfile SDF_OUTFILE
                            sdf output file
      -c CORES, --cores CORES
                            number of cores to use
      -p RUN_PLIP, --run-plip RUN_PLIP
                            Run PLIP?
      --victor VICTOR       Which victor to use: Victor, OpenVictor or Wictor
      -f IN_TABLE, --in-table IN_TABLE
                            CSV table input file, requires `name`, `smiles` and space-separated-`hit_names`

### Pipeline

Example:

```bash
fragmenstein pipeline \
                      --template template.min.pdb \
                      --input hits.sdf \
                      --n_cores $SLURM_CPUS_ON_NODE \
                      --suffix _frag \
                      --max_tasks 5000 \
                      --sw_databases REAL-Database-22Q1.smi.anon \
                      --combination_size 2 \
                      --workfolder /tmp/frag \
                      --timeout 600;
```

Usage:

    usage: fragmenstein pipeline [-h] -t TEMPLATE -i INPUT [-o OUTPUT] [-r RANKING] [-c CUTOFF] [-q QUICK] [-d SW_DIST] [-l SW_LENGTH] [-b SW_DATABASES [SW_DATABASES ...]] [-s SUFFIX]
                             [--workfolder WORKFOLDER] [--victor VICTOR] [-n N_CORES] [-m COMBINATION_SIZE] [-k TOP_MERGERS] [-e TIMEOUT] [-x MAX_TASKS] [-z BLACKLIST] [-j WEIGHTS] [-v]

    options:
      -h, --help            show this help message and exit
      -t TEMPLATE, --template TEMPLATE
                            Template PDB file
      -i INPUT, --input INPUT
                            Hits SDF file
      -o OUTPUT, --output OUTPUT
                            Output folder
      -r RANKING, --ranking RANKING
                            Ranking method
      -c CUTOFF, --cutoff CUTOFF
                            Joining cutoff
      -q QUICK, --quick QUICK
                            Quick reanimation
      -d SW_DIST, --sw_dist SW_DIST
                            SmallWorld distance
      -l SW_LENGTH, --sw_length SW_LENGTH
                            SmallWorld length
      -b SW_DATABASES [SW_DATABASES ...], --sw_databases SW_DATABASES [SW_DATABASES ...]
                            SmallWorld databases. Accepts multiple
      -s SUFFIX, --suffix SUFFIX
                            Suffix for output files
      --workfolder WORKFOLDER
                            Location to put the temp files
      --victor VICTOR       Which victor to use: Victor, OpenVictor or Wictor
      -n N_CORES, --n_cores N_CORES
                            Number of cores
      -m COMBINATION_SIZE, --combination_size COMBINATION_SIZE
                            Number of hits to combine in one step
      -k TOP_MERGERS, --top_mergers TOP_MERGERS
                            Max number of mergers to followup up on
      -e TIMEOUT, --timeout TIMEOUT
                            Timeout for each merger
      -x MAX_TASKS, --max_tasks MAX_TASKS
                            Max number of combinations to try in a batch
      -z BLACKLIST, --blacklist BLACKLIST
                            Blacklist file
      -j WEIGHTS, --weights WEIGHTS
                            JSON weights file
      -v, --verbose         verbose




* place the reference hits against themselves and gets the PLIP interactions
* combines the hits in given combination size, while skipping blacklisted named compounds.
* searches in [SmallWorld](sw.docking.org) the top N mergers
* places them and
* ranks them based on a customisable multiobjective function, which takes into account the PLIP interactions
     along with number of novel atoms (increase in risk & novelty).
 
This in effect reflects the pipeline I commonly use.

![pipeline](images/pipeline-01.png)

The values for the pipeline command are:

* `template`: The template, a polished PDB. The template must not contain a ligand in the site of interest,
                as Fragmenstein accepts other ligands (e.g. metals, cofactors etc.)
                and it is best to use a PyRosetta minimised template (i.e. one that has been through the ringer already).
* `hits`: The hits in sdf format. These need to have unique names.
* `output`: The output folder
* `suffix`: The suffix for the output files. Note that due to `max_tasks` there will be multiple sequential files for some steps.
* `quick`: Does not reattempt "reanimation" if it failed as the constraints are relaxed more and more the more deviation happens.
* `blacklist`: A file with a lines for each molecule name to not perform (say `hitA–hitZ`)
* `cutoff`: The joining cutoff in Ångström after which linkages will not be attempted (default is 5Å)
* `sw_databases`: See SmallWold or the [SmallWorld API in Python](https://github.com/matteoferla/Python_SmallWorld_API)
    for what datasets are available (e.g. 'Enamine-BB-Stock-Mar2022.smi.anon').
* `sw_length`: How many analogues for each query to keep
* `sw_dist`: The distance cutoff for the SmallWorld search
* `max_tasks`: To avoid memory issues, the pipeline performs a number of tasks (controlled via `max_tasks`)
    before processing them, to disable this use `--max_tasks 0`.
* `weights`: This is a JSON file that controls the ranking

This will minimise the first chain only stripping waters and heterogens

## Environment variables

To avoid having too many arguments, some default values can be set via environment variables.

A yaml file can be set in `$FRAGMENSTEIN_SETTINGS`,

The defaults can be seen in [settings.py](https://github.com/matteoferla/Fragmenstein/blob/9a026434dd27275c1b1f13b102b2e31e4c2b59ae/fragmenstein/settings.py)

In the yaml the values are lowercase, while as environment variables they are uppercase prefixed with `FRAGMENSTEIN_`.

Confusingly, for legacy reasons, the Victor argument, `.monster_throw_on_discard` controls `Monster.throw_on_discard`,
So `$FRAGMENSTEIN_MONSTER_THROW_ON_DISCARD` will not have an effect on `fragmenstein monster` run.

## Example slurm script

```bash
#!/bin/bash

# Some note for myself on usage
# submit via:
# export EXPERIMENT='_foo'; export TEMPLATE='protein.pdb'; export HITS='hits.sdf;
# export COMBO=2; export SLACK_WEBHOOK='';
# sbatch /mtn/someshared-drive/frag.slurm.sh;

#SBATCH --job-name=fragmenstein
#SBATCH --chdir=/mtn/someshared-drive/mferla
#SBATCH --output=/mtn/someshared-drive/logs/slurm-error_%x_%j.log
#SBATCH --error=/mtn/someshared-drive/logs/slurm-error_%x_%j.log
#SBATCH --clusters=clustername-this-is-cluster-specific
#SBATCH --partition=this-is-kind-of-the-group-and-is-cluster-specific
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --export=NONE,SLACK_WEBHOOK,COMBO,EXPERIMENT,TEMPLATE,HITS

# -------------------------------------------------------
# Some slurm fluff for logs

export SUBMITTER_HOST=$HOST
export HOST=$( hostname )
export USER=${USER:-$(users)}
export DATA=/mtn/someshared-drive
export HOME=$DATA/$USER; # homeless
source /etc/os-release;

echo "Running $SLURM_JOB_NAME ($SLURM_JOB_ID) as $USER in $HOST which runs $PRETTY_NAME submitted from $SUBMITTER_HOST"
echo "Request had cpus=$SLURM_JOB_CPUS_PER_NODE mem=$SLURM_MEM_PER_NODE tasks=$SLURM_NTASKS jobID=$SLURM_JOB_ID partition=$SLURM_JOB_PARTITION jobName=$SLURM_JOB_NAME"
echo "Started at $SLURM_JOB_START_TIME"
echo "job_pid=$SLURM_TASK_PID job_gid=$SLURM_JOB_GID topology_addr=$SLURM_TOPOLOGY_ADDR home=$HOME cwd=$PWD"

# -------------------------------------------------------
# Some conda fluff
# assuming you have an appropriate conda installed in `$HOME/conda-slurm`
source $HOME/conda-slurm/etc/profile.d/conda.sh
conda activate base;  # or whatever

pwd;
export EXPERIMENT=${EXPERIMENT:-'fragduo'}
export COMBO=${COMBO:-2}
export TEMPLATE=${TEMPLATE:-'template.pdb'}
export HITS=${HITS:-'hits.sdf'}
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX;
#export N_CORES=$(cat /proc/cpuinfo | grep processor | wc -l);
export N_CORES=$SLURM_CPUS_ON_NODE;

nice -19 fragmenstein pipeline \
                      --template $TEMPLATE \
                      --input $HITS \
                      --n_cores $(($N_CORES - 1)) \
                      --suffix _$EXPERIMENT \
                      --max_tasks 5000 \
                      --sw_databases REAL-Database-22Q1.smi.anon \
                      --combination_size $COMBO \
                      --workfolder /tmp/$EXPERIMENT \
                      --timeout 600;

rm -rf /tmp/$EXPERIMENT;

curl -X POST -H 'Content-type: application/json' --data '{"text":"'$EXPERIMENT' x2 complete"}' $SLACK_WEBHOOK
# -------------------------------------------------------

echo 'complete'

exit 0

```

