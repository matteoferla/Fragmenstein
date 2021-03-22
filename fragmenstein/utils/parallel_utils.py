import os
import logging

import dask
from dask.distributed import Client, LocalCluster

from fragmenstein.utils.config_manager import ConfigManager

# DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT="20" #Default is 10, try to set up

journal = logging.getLogger('Dask_Parallel')
journal.setLevel(logging.DEBUG)

DASK_CLIENT = None

def get_parallel_client(threads_per_worker=None, n_workers=None, memory_limit=None):
    global DASK_CLIENT
    if DASK_CLIENT is None:
        if n_workers is None:
            n_workers = ConfigManager.N_CPUS
            threads_per_worker = 1
        if memory_limit is None:
            memory_limit= ConfigManager.DASK_WORKER_MEMORY
            if memory_limit == "-1":
                from psutil import virtual_memory
                mem = virtual_memory()
                if mem.total is None:
                    raise ValueError("Error, memory was not determined")
                memory_limit="%dGB"%( (0.9 *mem.total/n_workers) // 2 ** 30)
        dask.config.set({'temporary_directory': os.path.join(ConfigManager.TMP_DIR, "dask")})
        if n_workers>1 or threads_per_worker>1:
            if ConfigManager.DISABLE_DASHBOARD:
                kwargs = {"dashboard_address":None}
            else:
                kwargs = {}
            DASK_CLIENT = Client( LocalCluster(threads_per_worker=threads_per_worker, n_workers=n_workers, memory_limit=memory_limit, **kwargs)) # dashboard_address=8787
        else:
            DASK_CLIENT = Client( LocalCluster(threads_per_worker=1, n_workers=1) )
        print(DASK_CLIENT, flush=True)
        if not ConfigManager.DISABLE_DASHBOARD:
            print(DASK_CLIENT.scheduler_info()['services'])
    return DASK_CLIENT


if __name__ == "__main__":
    print(get_parallel_client())