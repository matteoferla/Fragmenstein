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

    def __init__(self, random_seed=None, **victor_init_args):
        self.random_seed = random_seed
        self.victor_init_args = victor_init_args

        self.placed_victors = []

    def place(self,  smiles: str, number_runs: int= 10,
              long_name: str = 'ligand',
              merging_mode='expansion',
              atomnames: Optional[Dict[int, str]] = None,
              custom_map: Optional[Dict[str, Dict[int, int]]] = None,
              extra_ligand_constraint: Union[str] = None):
        """
        Places a followup (smiles) into the protein based upon the hits. Obtains number_runs solutions.
        :param smiles: smiles of followup, optionally covalent (_e.g._ ``*CC(=O)CCC``)
        :param number_runs: the number of Victor objects to create and place the smiles
        :param long_name: gets used for filenames so will get corrected
        :param merging_mode:
        :param atomnames: an optional dictionary that gets used by ``Params.from_smiles``
        :param custom_map: see Monster.place
        :param extra_ligand_constraint:
        :return:
        """
        assert number_runs >= 1, "Error, at least one placement required"
        current_randomState = None
        if self.random_seed:
            current_randomState = random.getstate()
            random.seed(self.random_seed)
        rseeds = [random.randint(0, int(1e6)) for i in range(number_runs)]
        if current_randomState is not None:
            random.setstate(current_randomState)

        i = 0
        victor = Victor(monster_random_seed=rseeds[i], **self.victor_init_args)
        victor.place(smiles, long_name=long_name+str(i), merging_mode= merging_mode,
                     custom_map=custom_map, atomnames=atomnames,
                     extra_ligand_constraint=extra_ligand_constraint)
        victor.runNumber = i
        ddG = victor.ddG
        monster = victor.monster
        self.placed_victors = [(ddG, victor)]

        for i in range(1, number_runs):
            monster.sample_new_conformation(rseeds[i])
            victor = Victor(monster_random_seed=rseeds[i], **self.victor_init_args)
            victor._prepare_args_for_placement(smiles=smiles, long_name = long_name+str(i), merging_mode= merging_mode,
                                               custom_map=custom_map,
                                               atomnames=atomnames, extra_ligand_constraint=extra_ligand_constraint)
            victor.monster = monster
            victor._calculate_placement_thermo()
            self.placed_victors += [(victor.ddG, victor)]
        self.placed_victors = sorted(self.placed_victors)

    def retrieve_best_victor(self):
        return self.placed_victors[0][1]

    def retrieve_scores(self):
        return list(zip(*self.placed_victors))[0]
