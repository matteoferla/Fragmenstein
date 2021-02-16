import os, sys
import tempfile
import warnings
from typing import Union, List

import numpy as np
from .rdkit_constrained_conformer_generation import gen_confs_constrained

from fragmenstein.external import ExternalToolImporter

warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore',category=FutureWarning)

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt


from collections import defaultdict


(DeLinker_test, frag_utils, rdkit_conf_parallel,
                data_prepare_data, example_utils)=  ExternalToolImporter.import_tool("DeLinker",
                                                      ["DeLinker_test", "frag_utils", "rdkit_conf_parallel", "data.prepare_data", "example_utils"])

DELINKER_ROOT = os.path.join(ExternalToolImporter.get_rootdir("DeLinker"))
DEEP_MODEL_FNAME = os.path.join(DELINKER_ROOT, "models/pretrained_DeLinker_model.pickle")
PAINS_DATA_FNAME = os.path.join(DELINKER_ROOT, "analysis/wehi_pains.csv")

class DeLinkerWrapper():

    DUMMY_SYMBOL = "*"

    def __init__(self, n_cores=6, gpu_id=0, interactive=False, random_seed=None):

        self.n_cores= n_cores
        self.gpu_id= gpu_id
        self.interactive= interactive
        self.random_seed = -1 if not random_seed else random_seed

    def load_mol(self, fnameOrMol, removeDummy=True, removeHs= True):
        #TODO: check if removeDummy interferes with covalent worheads

        if isinstance(fnameOrMol, str):
            mol = Chem.MolFromMolFile(fnameOrMol)
        else:
            mol = fnameOrMol
        if removeDummy:
          edit_mol = Chem.RWMol(mol)
          for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetSymbol() == DeLinkerWrapper.DUMMY_SYMBOL:
              edit_mol.RemoveAtom(i)

            mol= edit_mol.GetMol()
        if removeHs:
          mol= Chem.RemoveHs( mol )
        return mol


    def pick_closest_atoms(self, mol1, mol2, atom_idx1=None, atom_idx2=None, n_to_retrieve=None, cumm_idxs=True):
        '''

        :param mol1:
        :param mol2:
        :param atom_idx1:
        :param atom_idx2:
        :param n_to_retrieve:
        :param cumm_idxs: If true, indices of second molecule will start at N_ATOMS_1 position
        :return:
        '''

        if atom_idx1 and atom_idx2:
            return atom_idx1, atom_idx2

        # for c in mol1.GetConformers():
        #     print( c.GetPositions( ))

        combo = Chem.CombineMols(mol1, mol2)
        distance_matrix = Chem.Get3DDistanceMatrix(combo)
        n_atoms1 = mol1.GetNumAtoms()
        n_atoms2 = mol2.GetNumAtoms()

        if self.interactive:
            mol_with_idxs = example_utils.mol_with_atom_index(combo_2d)
            plt.imshow(Draw.MolsToGridImage([mol_with_idxs], molsPerRow=1))
            plt.show()

        distance_matrix = distance_matrix[:n_atoms1, n_atoms2:].copy() #TODO: Esto tiene que estar mal

        if atom_idx1:
            distance_matrix_tmp = np.ones_like(distance_matrix) * sys.maxsize
            distance_matrix_tmp[atom_idx1, :] = distance_matrix
            distance_matrix = distance_matrix_tmp

        if atom_idx2:
            distance_matrix_tmp = np.ones_like(distance_matrix) * sys.maxsize
            distance_matrix_tmp[:, atom_idx1] = distance_matrix
            distance_matrix = distance_matrix_tmp

        unravel_idxs= np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)
        selected_pairs=[]
        k = 0
        if not n_to_retrieve:
            n_to_retrieve = len(unravel_idxs[0])
        for i,j in zip(* unravel_idxs ):
            k += 1
            if k > n_to_retrieve:
                break
            selected_pairs.append( (int(i), int(j+ int(cumm_idxs)*n_atoms1)) )

        return selected_pairs

    def _link_molecule_pairs_given_atom_idxs(self, frag_1_mol: Chem.Mol, frag_2_mol: Chem.Mol, atom_idx_1:int, atom_idx_2:int) -> List[Chem.Mol]:
        '''
        :param frag_1_mol: A molecule to be linked to frag_2_mol
        :param frag_2_mol: A molecule to be linked to frag_1_mol
        :param atom_idx_1:
        :param atom_idx_2:
        :return:
        '''

        combo_no_exit = Chem.CombineMols(frag_1_mol, frag_2_mol)
        combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))  #TODO: Check if *.* can be used due to DUMMY symbol

        combo_2d = Chem.Mol(combo)
        AllChem.Compute2DCoords(combo_2d)
        num_heavy_atoms = combo.GetNumHeavyAtoms()

        edcombo = Chem.EditableMol(combo)

        num_heavy_atoms = combo.GetNumHeavyAtoms()

        edcombo.AddBond(num_heavy_atoms, atom_idx_1, order=Chem.rdchem.BondType.SINGLE)
        edcombo.AddBond(num_heavy_atoms + 1, atom_idx_2, order=Chem.rdchem.BondType.SINGLE)
        editedcombo = edcombo.GetMol()

        AllChem.Compute2DCoords(editedcombo)

        AllChem.Compute2DCoords(editedcombo)
        Chem.SanitizeMol(editedcombo)

        if self.interactive:
            plt.imshow( Draw.MolsToGridImage([editedcombo], molsPerRow=1) ); plt.show()


        mol_to_link = edcombo.GetMol()
        Chem.SanitizeMol(mol_to_link)

        # Convert exit vectors to carbons for conformer generation
        du = Chem.MolFromSmiles('*')
        mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link, du, Chem.MolFromSmiles('C'), True)[0]
        Chem.SanitizeMol(mol_to_link_carbon)
        # Generate conformer
        mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
        AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit, randomseed=self.random_seed)
        mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)

        # Add this conformer to the two unlinked fragments
        conf = mol_to_link.GetConformer()
        ref_conf = mol_to_link_carbon.GetConformer()
        for i in range(mol_to_link_carbon.GetNumAtoms()):
            pos = list(ref_conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, pos)
        conf.SetId(0)
        _ = mol_to_link.AddConformer(conf)

        if self.interactive:
            plt.imshow( Draw.MolsToGridImage([mol_to_link], molsPerRow=1) ); plt.show()


        # Get distance and angle between fragments
        dist, ang = frag_utils.compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))

        print(Chem.MolToSmiles(mol_to_link), dist, ang)

        cwdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            data_path = "./fragments_test_data.txt"
            with open(data_path, 'w') as f:
                f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
            raw_data = data_prepare_data.read_file(data_path)
            data_prepare_data.preprocess(raw_data, "zinc", "fragments_test", True)

            if not self.gpu_id:
                os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
            else:
                os.environ['CUDA_VISIBLE_DEVICES'] = str(self.gpu_id)

            # Arguments for DeLinker_to_remove
            args = defaultdict(None)
            args['--dataset'] = 'zinc'
            args['--config'] = '{"generation": true, \
                                               "batch_size": 1, \
                                               "number_of_generation_per_valid": 50, \
                                               "min_atoms": 2, "max_atoms": 6, \
                                               "train_file": "molecules_fragments_test.json", \
                                               "valid_file": "molecules_fragments_test.json", \
                                               "output_name": "DeLinker_example_generation.smi"}'
            args['--freeze-graph-model'] = False
            args['--restore'] = DEEP_MODEL_FNAME

            # Setup model and generate molecules
            model = DeLinker_test.DenseGGNNChemModel(args)
            generated_smiles = model.generate_new_graphs(model.valid_data, return_list=True)
            del model
            # print(generated_molecules)
            # input("enter")
            # model.train()
            # # Free up some memory
            # del model
            #
            # # Load molecules
            # generated_smiles = frag_utils.read_triples_file("./DeLinker_example_generation.smi")
        os.chdir(cwdir)
        in_mols = [smi[1] for smi in generated_smiles]
        frag_mols = [smi[0] for smi in generated_smiles]
        gen_mols = [smi[2] for smi in generated_smiles]

        du = Chem.MolFromSmiles('*')
        clean_frags = [Chem.MolToSmiles(Chem.RemoveHs(
            AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi), du, Chem.MolFromSmiles('[H]'), True)[0])) for smi in
                       frag_mols]

        print("Deep learning done!")

        # Check valid
        results = []
        for in_mol, frag_mol, gen_mol, clean_frag in zip(in_mols, frag_mols, gen_mols, clean_frags):
            if len(Chem.MolFromSmiles(gen_mol).GetSubstructMatch(Chem.MolFromSmiles(clean_frag))) > 0:
                results.append([in_mol, frag_mol, gen_mol, clean_frag])

        print("Number of generated SMILES: \t%d" % len(generated_smiles))
        print("Number of valid SMILES: \t%d" % len(results))
        print("%% Valid: \t\t\t%.2f%%" % (len(results) / len(generated_smiles) * 100))

        import re

        linkers = list(map(frag_utils.get_linker, [Chem.MolFromSmiles(m[2]) for m in results],
                           [Chem.MolFromSmiles(m[3]) for m in results], [m[1] for m in results]
                           ))
        # Standardise linkers
        for i, linker in enumerate(linkers):
            if linker == "":
                continue
            linker = Chem.MolFromSmiles(re.sub('[0-9]+\*', '*', linker))
            Chem.rdmolops.RemoveStereochemistry(linker)
            linkers[i] = Chem.MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker))
        # Update results
        for i in range(len(results)):
            results[i].append(linkers[i])

        print("Done")

        # Create dictionary of results
        results_dict = {}
        for res in results:
            if res[0] + '.' + res[
                1] in results_dict:  # Unique identifier - starting fragments and original molecule
                results_dict[res[0] + '.' + res[1]].append(tuple(res))
            else:
                results_dict[res[0] + '.' + res[1]] = [tuple(res)]
        # Check uniqueness
        print("Unique molecules: %.2f%%" % (frag_utils.unique(results_dict.values()) * 100))

        # Check if molecules pass 2D filters
        filters_2d = frag_utils.calc_filters_2d_dataset(results,
                                                        pains_smarts_loc=PAINS_DATA_FNAME,
                                                        n_cores= self.n_cores)

        results_filt = []
        for res, filt in zip(results, filters_2d):
            if filt[0] and filt[1] and filt[2]:
                results_filt.append(res)

        print("Pass all 2D filters: \t\t\t\t%.2f%%" % (len(results_filt) / len(results) * 100))
        print("Valid and pass all 2D filters: \t\t\t%.2f%%" % (len(results_filt) / len(generated_smiles) * 100))
        print(
            "Pass synthetic accessibility (SA) filter: \t%.2f%%" % (
                len([f for f in filters_2d if f[0]]) / len(filters_2d) * 100))
        print("Pass ring aromaticity filter: \t\t\t%.2f%%" % (
            len([f for f in filters_2d if f[1]]) / len(filters_2d) * 100))
        print(
            "Pass SA and ring filters: \t\t\t%.2f%%" % (
                len([f for f in filters_2d if f[0] and f[1]]) / len(filters_2d) * 100))
        print("Pass PAINS filters: \t\t\t\t%.2f%%" % (len([f for f in filters_2d if f[2]]) / len(filters_2d) * 100))

        # Get unique molecules
        print("Number molecules passing 2D filters:\t\t%d" % len(results_filt))
        results_filt_unique = example_utils.unique_mols(results_filt)
        print("Number unique molecules passing 2D filters:\t%d" % len(results_filt_unique))

        new_gen_mols = [x[2] for x in results_filt_unique]

        args = dict(smiles_full_mols=[new_gen_mols], sdffile=None, reference_confs=[mol_to_link], maxconfs=50,
                    rms_threshold=0.35, energy=10,
                    tdist=0.75, smi_frags=[frag_mols], numcores=self.n_cores, jpsettings=False, verbose=True)


        confs = gen_confs_constrained(**args)

        print( confs[0], type(confs[0]) )
        return confs

    def link_molecule_pair(self, fnameOrMol_1: Union[str,Chem.Mol], fnameOrMol_2:Union[str, Chem.Mol], n_attempts=3) -> List[Chem.Mol]:
        '''
        :param fnameOrMol_1:
        :param fnameOrMol_2:
        :param atom_idx_1:
        :param atom_idx_2:
        :return:
        '''
        frag_1_mol = self.load_mol(fnameOrMol_1)

        frag_2_mol = self.load_mol(fnameOrMol_2)

        linked_molecules = []
        for atom_idx_1, atom_idx_2 in self.pick_closest_atoms(frag_1_mol, frag_2_mol ):
            linked_molecules += self._link_molecule_pairs_given_atom_idxs( frag_1_mol, frag_2_mol, atom_idx_1, atom_idx_2)
            n_attempts-=1
            if n_attempts<1:
                break

        return linked_molecules

# How many cores for multiprocessing
n_cores = 6
# Whether to use GPU for generating molecules with DeLinker_to_remove


frag_1_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x1458.mol"
frag_2_path="/home/ruben/oxford/tools/Fragmenstein/fragmenstein/mpro/data/hit_mols/Mpro-x0104.mol"

dlw= DeLinkerWrapper()
sp= dlw.link_molecule_pair(dlw.load_mol(frag_1_path), dlw.load_mol(frag_2_path))
print( sp )
# sys.exit(0)


# frag_1_sdf = dlw.load_mol(frag_1_path)
# frag_1_smi = Chem.MolToSmiles(frag_1_sdf)
#
# frag_2_sdf = dlw.load_mol(frag_2_path)
# frag_2_smi = Chem.MolToSmiles(frag_2_sdf)
#
# # print( Chem.MolToPDBBlock( frag_1_sdf ) )
# # img = Draw.MolsToGridImage([Chem.MolFromSmiles(frag_1_smi), Chem.MolFromSmiles(frag_2_smi)], molsPerRow=2, subImgSize=(300, 300))
# # plt.imshow(img); plt.show()
#
# combo_no_exit = Chem.CombineMols(frag_1_sdf, frag_2_sdf)
# combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))
#
# # for i in range(combo.GetNumAtoms()):
# #   atom= combo.GetAtomWithIdx(i)
# #   print(i, atom.GetSymbol(), atom.GetExplicitValence() )
#
# combo_2d = Chem.Mol(combo)
# AllChem.Compute2DCoords(combo_2d)
# mol_with_idxs= example_utils.mol_with_atom_index(combo_2d)
# plt.imshow( Draw.MolsToGridImage([mol_with_idxs], molsPerRow=1) ); plt.show()
#
#
# atom_idx_1 = 26
# atom_idx_2 = 12
#
#
# edcombo = Chem.EditableMol(combo)
#
# num_heavy_atoms = combo.GetNumHeavyAtoms()
#
# edcombo.AddBond(num_heavy_atoms, atom_idx_1, order=Chem.rdchem.BondType.SINGLE)
# edcombo.AddBond(num_heavy_atoms+1, atom_idx_2, order=Chem.rdchem.BondType.SINGLE)
# editedcombo = edcombo.GetMol()
#
#
# AllChem.Compute2DCoords(editedcombo)
#
#
#
# AllChem.Compute2DCoords(editedcombo)
# Chem.SanitizeMol(editedcombo )
#
# # plt.imshow( Draw.MolsToGridImage([editedcombo], molsPerRow=1) ); plt.show()
#
#
# mol_to_link = edcombo.GetMol()
# Chem.SanitizeMol(mol_to_link)
#
# # Convert exit vectors to carbons for conformer generation
# du = Chem.MolFromSmiles('*')
# mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link,du,Chem.MolFromSmiles('C'),True)[0]
# Chem.SanitizeMol(mol_to_link_carbon)
# # Generate conformer
# mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
# AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit, randomseed=42)
# mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)
#
# # Add this conformer to the two unlinked fragments
# conf = mol_to_link.GetConformer()
# ref_conf = mol_to_link_carbon.GetConformer()
# for i in range(mol_to_link_carbon.GetNumAtoms()):
#     pos = list(ref_conf.GetAtomPosition(i))
#     conf.SetAtomPosition(i, pos)
# conf.SetId(0)
# _ = mol_to_link.AddConformer(conf)
#
# # plt.imshow( Draw.MolsToGridImage([mol_to_link], molsPerRow=1) ); plt.show()
#
#
# # Get distance and angle between fragments
# dist, ang = frag_utils.compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))
#
# print(Chem.MolToSmiles(mol_to_link), dist, ang )
#
# cwdir= os.getcwd()
# with tempfile.TemporaryDirectory() as tmpdirname:
#   os.chdir(tmpdirname)
#   data_path = "./fragments_test_data.txt"
#   with open(data_path, 'w') as f:
#     f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
#   raw_data = data_prepare_data.read_file(data_path)
#   data_prepare_data.preprocess(raw_data, "zinc", "fragments_test", True)
#
#   if not use_gpu:
#       os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
#
#   # Arguments for DeLinker_to_remove
#   args = defaultdict(None)
#   args['--dataset'] = 'zinc'
#   args['--config'] = '{"generation": true, \
#                        "batch_size": 1, \
#                        "number_of_generation_per_valid": 50, \
#                        "min_atoms": 2, "max_atoms": 6, \
#                        "train_file": "molecules_fragments_test.json", \
#                        "valid_file": "molecules_fragments_test.json", \
#                        "output_name": "DeLinker_example_generation.smi"}'
#   args['--freeze-graph-model'] = False
#   args['--restore'] = DEEP_MODEL_FNAME
#
#   # Setup model and generate molecules
#   model = DeLinker_test.DenseGGNNChemModel(args)
#   model.train()
#   # Free up some memory
#   del model
#
#   # Load molecules
#   generated_smiles = frag_utils.read_triples_file("./DeLinker_example_generation.smi")
# os.chdir(cwdir)
# in_mols = [smi[1] for smi in generated_smiles]
# frag_mols = [smi[0] for smi in generated_smiles]
# gen_mols = [smi[2] for smi in generated_smiles]
#
# du = Chem.MolFromSmiles('*')
# clean_frags = [Chem.MolToSmiles(Chem.RemoveHs(AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi),du,Chem.MolFromSmiles('[H]'),True)[0])) for smi in frag_mols]
#
# print("Deep learning done!")
#
#
# # Check valid
# results = []
# for in_mol, frag_mol, gen_mol, clean_frag in zip(in_mols, frag_mols, gen_mols, clean_frags):
#     if len(Chem.MolFromSmiles(gen_mol).GetSubstructMatch(Chem.MolFromSmiles(clean_frag)))>0:
#         results.append([in_mol, frag_mol, gen_mol, clean_frag])
#
# print("Number of generated SMILES: \t%d" % len(generated_smiles))
# print("Number of valid SMILES: \t%d" % len(results))
# print("%% Valid: \t\t\t%.2f%%" % (len(results)/len(generated_smiles)*100))
#
# import re
#
# linkers = list(map(frag_utils.get_linker, [Chem.MolFromSmiles(m[2]) for m in results],
#                    [Chem.MolFromSmiles(m[3]) for m in results], [m[1] for m in results]
#                    ))
# # Standardise linkers
# for i, linker in enumerate(linkers):
#   if linker == "":
#     continue
#   linker = Chem.MolFromSmiles(re.sub('[0-9]+\*', '*', linker))
#   Chem.rdmolops.RemoveStereochemistry(linker)
#   linkers[i] = Chem.MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker))
# # Update results
# for i in range(len(results)):
#   results[i].append(linkers[i])
#
# print("Done")
#
# # Create dictionary of results
# results_dict = {}
# for res in results:
#     if res[0]+'.'+res[1] in results_dict: # Unique identifier - starting fragments and original molecule
#         results_dict[res[0]+'.'+res[1]].append(tuple(res))
#     else:
#         results_dict[res[0]+'.'+res[1]] = [tuple(res)]
# # Check uniqueness
# print("Unique molecules: %.2f%%" % (frag_utils.unique(results_dict.values())*100))
#
# # Check if molecules pass 2D filters
# filters_2d = frag_utils.calc_filters_2d_dataset(results,
#                                                 pains_smarts_loc= PAINS_DATA_FNAME,
#                                                 n_cores=n_cores)
#
# results_filt = []
# for res, filt in zip(results, filters_2d):
#   if filt[0] and filt[1] and filt[2]:
#     results_filt.append(res)
#
# print("Pass all 2D filters: \t\t\t\t%.2f%%" % (len(results_filt) / len(results) * 100))
# print("Valid and pass all 2D filters: \t\t\t%.2f%%" % (len(results_filt) / len(generated_smiles) * 100))
# print(
#   "Pass synthetic accessibility (SA) filter: \t%.2f%%" % (len([f for f in filters_2d if f[0]]) / len(filters_2d) * 100))
# print("Pass ring aromaticity filter: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[1]]) / len(filters_2d) * 100))
# print(
#   "Pass SA and ring filters: \t\t\t%.2f%%" % (len([f for f in filters_2d if f[0] and f[1]]) / len(filters_2d) * 100))
# print("Pass PAINS filters: \t\t\t\t%.2f%%" % (len([f for f in filters_2d if f[2]]) / len(filters_2d) * 100))
#
# # Get unique molecules
# print("Number molecules passing 2D filters:\t\t%d" % len(results_filt))
# results_filt_unique = example_utils.unique_mols(results_filt)
# print("Number unique molecules passing 2D filters:\t%d" % len(results_filt_unique))
#
#
# new_gen_mols= [ x[2] for x in results_filt_unique ]
# args= dict(smiles_full_mols= [new_gen_mols] , sdffile= None, reference_confs= [mol_to_link], maxconfs=50, rms_threshold=0.35, energy=10,
#           tdist=0.75, smi_frags= [frag_mols], numcores=n_cores, jpsettings=False, verbose=True )
#
# from joblib import dump
# dumpFname= os.path.join(THIS_SCRIPT_PARENTDIR, "args.pkl")
# print( dumpFname )
# dump(args, dumpFname)
#
# confs = gen_confs_constrained ( [gen_mols] , None, [mol_to_link], maxconfs=10, rms_threshold=0.35, energy=10,
#           tdist=0.75, smi_frags= [frag_mols], numcores=n_cores, jpsettings=False)
#
#
# print(confs)
#
# dumpFname= os.path.join(THIS_SCRIPT_PARENTDIR, "molecules_generated.pkl")
# print( dumpFname )
# dump(confs, dumpFname)