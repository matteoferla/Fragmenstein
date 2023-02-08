import enum


class VictorMinMode(enum.Enum):
    """Victor mode for minimisation
    TODO this is not used yet
    """
    NULL = -1  #: no minimisation at all
    IGOR = 0  #: use Igor's minimisation (PyRosetta)
    IGOR_NODISK = 1  #: use Igor's minimisation without writing to disk
    FRITZ = 2  #: use Fritz's minimisation (openMM)
