import copy
import random
from typing import Optional, Dict, Union

from rdkit import Chem

from fragmenstein import Victor


class MultiVictorPlacement():
    """
    MultiVictorPlacement is a class to execute multiple times victor.place over the same smile using different seeds.
    The idea is to be able to do something more like docking for those cases in which the inspirational hits do not
    explain an important part of the molecule
    """

    def __init__(self, randomSeed=None, **victor_init_args):
        self.randomSeed = randomSeed
        self.victor_init_args = victor_init_args

        self.placed_victors = []

    def place(self,  smiles: str, number_trials: int= 10,
              long_name: str = 'ligand',
              merging_mode='none_permissive',
              atomnames: Optional[Dict[int, str]] = None,
              extra_ligand_constraint: Union[str] = None):
        current_randomState = None
        if self.randomSeed:
            current_randomState = random.getstate()
            random.seed(self.randomSeed)
        rseeds = [random.randint(0, int(1e6)) for i in range(number_trials)]
        if current_randomState is not None:
            random.setstate(current_randomState)

        def placeFunction(i):
            victor = Victor(randomSeed=rseeds[i], **self.victor_init_args)
            victor.place(smiles=smiles, long_name = long_name+str(i), merging_mode= merging_mode, atomnames=atomnames,
                         extra_ligand_constraint=extra_ligand_constraint)
            victor.runNumber = i
            return victor.ddG, victor

        self.placed_victors = sorted((map(placeFunction, range(number_trials))))

    def retrieve_best_victor(self):
        return self.placed_victors[0][1]

    def retrieve_scores(self):
        return list(zip(*self.placed_victors))[0]