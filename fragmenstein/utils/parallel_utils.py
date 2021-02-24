import os
import logging

import dask
from distributed import Client

from fragmenstein.utils.config_manager import ConfigManager

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
            if memory_limit <0:
                from psutil import virtual_memory
                mem = virtual_memory()
                if mem.total is None:
                    print(mem)
                    raise ValueError("Error, memory was not determined")
                memory_limit="%dGB"%( (0.9 *mem.total/n_workers) // 2 ** 30)
        dask.config.set({'temporary_directory': os.path.join(ConfigManager.TMP_DIR, "dask")})
        if n_workers>1 or threads_per_worker>1:
            DASK_CLIENT = Client(threads_per_worker=threads_per_worker, n_workers=n_workers, memory_limit=memory_limit)
        else:
            DASK_CLIENT = Client(threads_per_worker=1, n_workers=1)

    return DASK_CLIENT


if __name__ == "__main__":
    print(get_parallel_client())