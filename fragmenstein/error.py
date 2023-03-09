from typing import Sequence, Union
from rdkit import Chem

class FragmensteinError(Exception):
    pass

class ShoddyCodeError(FragmensteinError):
    """
    This is a placeholder for an error that should not happen!
    """
    pass

class FullOverlapError(FragmensteinError):
    def __init__(self,
                 message='Full overlap with no constructive outcome.',
                 hits: Sequence[Union[Chem.Mol, str]] = ()):
        self.message = message
        self.hits = hits
class DistanceError(FragmensteinError):
    """
    This was formerly ConnectionError for comedy.
    """
    def __init__(self, distance: float = float('nan'), hits: Sequence[Union[Chem.Mol, str]] = (), message=None):
        self.message = message
        self.distance = distance
        self.hits = hits
    def __str__(self):
        hit_names = '+'.join([h.GetProp('_Name') if isinstance(h, Chem.Mol) else h for h in self.hits])
        distance = f'{self.distance:.2f} Å' if str(self.distance) != 'nan' else ''
        parts = [self.message, hit_names, distance]
        return ' '.join([p for p in parts if p])

class RectificationError(FragmensteinError):

    def __init__(self, message=None, mol=None):
        self.message = message
        self.mol = mol
    def __str__(self):
        return self.message
