#  ####################
  #
  # Example 1
  # Simple condor job description file
  #
  ####################

Executable      = /data/xchem-fragalysis/sanchezg/app/miniconda3_2/envs/Fragmenstein/bin/python
Arguments       = -m examples.protocols_bricsFragmenstein --n_cpus 1 -d /data/xchem-fragalysis/sanchezg/oxford/myProjects/diamondCovid/data/nsp13/aligned -o /data/xchem-fragalysis/sanchezg/oxford/tools/Fragmenstein/output_Site_in -f x0041_0A x0176_0A x0276_0A
Universe        = vanilla
Output          = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).out
Error           = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).err
Log             = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).log

environment = "EXTERNAL_TOOLS_CONFIG_FILE=examples/external_config.json DASK_WORKER_MEMORY=4GB"

request_cpus   = 1
#request_memory = 4096
#request_disk   = 262144
#request_gpus    = 1

#requirements = (TARGET.Machine == "pulsar-exec-node-gpu-2.xchem.novalocal")

Queue
