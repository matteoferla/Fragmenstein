__version__ = "0.12.3"

from typing import Dict
def get_versions() -> Dict[str, str]:
    """
    Return a dict of versions of os, python, fragmenstein etc.
    """
    import pkg_resources, sys, platform

    get_version = lambda name: pkg_resources.get_distribution(name).version

    return dict(python=sys.version, os_type=platform.system(), arc=platform.machine(),
                fragmenstein=get_version("fragmenstein"),
                pyrosetta=get_version("pyrosetta"),
                rdkit=get_version("rdkit"),
                molecular_rectifier=get_version('molecular-rectifier'),
                rdkit_to_params=get_version('rdkit-to-params')
                )
