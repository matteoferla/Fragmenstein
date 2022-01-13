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

    def place(self,  smiles: str, number_trials: int= 10,
              long_name: str = 'ligand',
              merging_mode='none_permissive',
              atomnames: Optional[Dict[int, str]] = None,
              extra_ligand_constraint: Union[str] = None):
        assert number_trials >= 1, "Error, at least one placement required"
        current_randomState = None
        if self.random_seed:
            current_randomState = random.getstate()
            random.seed(self.random_seed)
        rseeds = [random.randint(0, int(1e6)) for i in range(number_trials)]
        if current_randomState is not None:
            random.setstate(current_randomState)

        def placeFunction(i):
            victor = Victor(random_seed=rseeds[i], **self.victor_init_args)
            victor.place(smiles=smiles, long_name = long_name+str(i), merging_mode= merging_mode, atomnames=atomnames,
                         extra_ligand_constraint=extra_ligand_constraint)
            victor.runNumber = i
            return victor.ddG, victor

        ddG, victor = placeFunction(0)
        monster = victor.monster
        self.placed_victors = [ddG, victor]

        monster.sample_new_conformation()

        i=1
        victor = Victor(random_seed=rseeds[i], **self.victor_init_args)
        victor._prepare_args_for_placement(smiles=smiles, long_name = long_name+str(i), merging_mode= merging_mode,
                                           atomnames=atomnames, extra_ligand_constraint=extra_ligand_constraint)
        victor.monster = monster
        victor._calculate_placement_minimizeMonster()
        self.placed_victors = [victor.ddG, victor]

        # self.placed_victors = map(placeFunction, range(number_trials))
        self.placed_victors = sorted(self.placed_victors)

    def retrieve_best_victor(self):
        return self.placed_victors[0][1]

    def retrieve_scores(self):
        return list(zip(*self.placed_victors))[0]

if __name__ == '__main__':
    import sys, os
    from rdkit import Chem
    import pyrosetta

    test_path="/home/sanchezg/oxford/tools/Fragmenstein/test_mols/"
    pyrosetta.init(
        extra_options='-no_optH false -mute all -ex1 -ex2 -ignore_unrecognized_res false -load_PDB_components false -ignore_waters false')
    to_place = Chem.MolFromMolFile(test_path+'placed_example1.mol')
    pdb_filename = test_path+'apo_example1.pdb'
    smiles = Chem.MolToSmiles(to_place)
    hits = [ Chem.MolFromMolFile(os.path.join(test_path, basename)) for basename in ["x0032_0A.mol", "x0103_0A.mol"]]
    mv = MultiVictorPlacement(hits=hits, pdb_filename=pdb_filename)
    mv.place(smiles)
    print(mv.retrieve_best_victor())

# class MultiVictorPlacement():
#     """
#     MultiVictorPlacement is a class to execute multiple times victor.place over the same smile using different seeds.
#     The idea is to be able to do something more like docking for those cases in which the inspirational hits do not
#     explain an important part of the molecule
#     """
#
#     def __init__(self, random_seed=None, **victor_init_args):
#         self.random_seed = random_seed
#         self.victor_init_args = victor_init_args
#
#         self.placed_victors = []
#
#     def place(self,  smiles: str, number_trials: int= 10,
#               long_name: str = 'ligand',
#               merging_mode='none_permissive',
#               atomnames: Optional[Dict[int, str]] = None,
#               extra_ligand_constraint: Union[str] = None):
#         current_randomState = None
#         if self.random_seed:
#             current_randomState = random.getstate()
#             random.seed(self.random_seed)
#         rseeds = [random.randint(0, int(1e6)) for i in range(number_trials)]
#         if current_randomState is not None:
#             random.setstate(current_randomState)
#
#         def placeFunction(i):
#             victor = Victor(random_seed=rseeds[i], **self.victor_init_args)
#             victor.place(smiles=smiles, long_name = long_name+str(i), merging_mode= merging_mode, atomnames=atomnames,
#                          extra_ligand_constraint=extra_ligand_constraint)
#             victor.runNumber = i
#             return victor.ddG, victor
#
#         self.placed_victors = sorted((map(placeFunction, range(number_trials))))
#
#     def retrieve_best_victor(self):
#         return self.placed_victors[0][1]
#
#     def retrieve_scores(self):
#         return list(zip(*self.placed_victors))[0]