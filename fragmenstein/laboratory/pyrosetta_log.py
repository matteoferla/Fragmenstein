## Copied from https://github.com/matteoferla/pyrosetta_scripts/tree/main/init_helper

import io
import logging
import pyrosetta
import re
from typing import Union, List


def configure_logger() -> logging.Logger:
    """
    The function `get_logger`, simply adds a stringIO handler to the log and captures the log,
    thus making it easier to use.
    The function `get_log_entries`, spits out entries of a given level.

    :return: logger
    """
    pyrosetta.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.INFO)  # default = logging.WARNING
    stringio = io.StringIO()
    handler = logging.StreamHandler(stringio)
    handler.setLevel(logging.INFO)
    # handler.set_name('stringio')
    handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
    logger.addHandler(handler)
    return logger


def get_log_entries(levelname: Union[str, int] = logging.INFO) -> List[str]:
    """
    Get a list of all entries in log at a given level.
    levelname can be either an int (``logging.INFO`` etc. are numbers multiples of 10 in increasing severity)
    or a string of the level.
    Note that it is very crude: if INFO is requested, ERROR is not shown!

    :param levelname: int for the level number or str of the name
    :return: List of str
    """
    if isinstance(levelname, int):
        # logging.INFO is actually an int, not an enum
        levelname = logging.getLevelName(levelname)
    stringio = logging.getLogger("rosetta").handlers[0].stream
    return re.findall(f'(\[.*\] {levelname} - [\w\W]*)', stringio.getvalue())
