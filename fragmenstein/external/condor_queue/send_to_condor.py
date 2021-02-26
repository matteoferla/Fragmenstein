import argparse
import os
import tempfile

from subprocess import check_call

DEFAULT_OGS_DIR = "/data/xchem-fragalysis/sanchezg/logs/"
SUBMITS_DIR = "/data/xchem-fragalysis/sanchezg/submits/"

BASH_TEMPLATE='''###################
%(env_vars)s
%(cmd)s

###################
'''

CONDOR_TEMPLATE='''###################
Executable      = /usr/bin/bash
Arguments       = %(bash_tmpFile)s
Universe        = vanilla
Output          = %(logdirs)s/$(Cluster).$(Process).out
Error           = %(logdirs)s/$(Cluster).$(Process).err
Log             = %(logdirs)s/$(Cluster).$(Process).log

request_cpus   = %(ncpus)s
%(additional_requirements)s

Queue

###################

'''


'''
#additional_requirements example:

#request_memory = 4096
#request_disk   = 262144
#request_gpus    = 1

#requirements = (TARGET.Machine == "pulsar-exec-node-gpu-2.xchem.novalocal")

###################
'''

parser = argparse.ArgumentParser("utility to send commands to condor queue")

parser.add_argument("--ncpus", type=int, required=True, help="number of cpus")
parser.add_argument("--memory", type=int, required=False, default=None, help="Total memory in MB. Default %(default)s")
parser.add_argument("--gpus", type=int, required=False, default=None, help="Number of GPUs")
# parser.add_argument("--bindir", type=str, required=False, default=None, help="directory where the binary lives") #TODO
parser.add_argument("--logdirs", type=str, required=False, default=DEFAULT_OGS_DIR, help="Logs directory")
parser.add_argument("--env_vars", type=str, nargs="+", required=False, default=[], help="enviramental variables")
parser.add_argument("--print", action="store_true", help="print files instead submitting")

parser.add_argument("cmd", type=str, help="commands between \"\"")

args = vars(parser.parse_args())
print( args )

if args:

    env_vars =""
    for env_var in args["env_vars"]:
        env_vars+= ("export " +env_var+ "\n")
    args["env_vars"]= env_vars
    bash_str = BASH_TEMPLATE%(args)
    additional_requirements=""
    if args["memory"]:
        additional_requirements += "request_memory = %d\n"%args["memory"]
    if args["gpus"]:
        additional_requirements += "request_gpus = %d\n"%args["gpus"]
    args["additional_requirements"] = additional_requirements
    if args["print"]:
        print("***********************************")
        print(bash_str )
        print("***********************************")

    with tempfile.TemporaryDirectory(dir=SUBMITS_DIR) as tmp:
        with tempfile.NamedTemporaryFile(dir = SUBMITS_DIR, suffix="_launch.sh", delete=False) as tmpfile:
            bash_tmpFile = tmpfile.name
            tmpfile.write( bash_str)
        args["bash_tmpFile"] = bash_tmpFile
        condor_str = CONDOR_TEMPLATE%(args)

        if args["print"]:
            print("***********************************")
            print(condor_str)
            print("***********************************")
        else:
            condor_tmpFile = os.path.join(tmp, "launch.condor")
            with open( condor_tmpFile, "w") as tmpfile:
                tmpfile.write( condor_str)

            check_call(["condor_submit", condor_tmpFile ])
'''
python -m fragmenstein.external.condor_queue.send_to_condor --ncpus 4 "python -c 'print(0)'"
'''