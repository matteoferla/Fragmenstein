from abc import ABC, abstractmethod


class LoadInput_base(ABC):


    def __init__(self):
        self._fragments_dict = None


    @abstractmethod
    def prepare_fragments(self):
        raise  NotImplementedError

    @property
    def fragments_dict(self):
        if not self._fragments_dict:
            self._fragments_dict = self.prepare_fragments()
        return self._fragments_dict

    @property
    def fragments(self):
        return iter( self.fragments_dict.values())


