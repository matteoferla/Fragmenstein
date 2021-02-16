  ####################
  #
  # Example 1
  # Simple condor job description file
  #
  ####################

Executable      = /usr/bin/bash
Arguments       = /data/xchem-fragalysis/sanchezg/oxford/run_simulation.sh
Universe        = vanilla
Output          = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).out.txt
Error           = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).err.txt
Log             = /data/xchem-fragalysis/sanchezg/logs/$(Cluster).$(Process).log.txt

#request_cpus   = 1
#request_memory = 4096
#request_disk   = 262144
request_gpus    = 1

requirements = (TARGET.Machine == "pulsar-exec-node-gpu-2.xchem.novalocal")

Queue




