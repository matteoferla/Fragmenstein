import os, sys, re
import tempfile
import warnings
from typing import Union, List

import numpy as np
from rdkit import RDLogger

from fragmenstein.external import ExternalToolImporter

warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore',category=FutureWarning)

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from matplotlib import pyplot as plt


from collections import defaultdict




DELINKER_ROOT = os.path.join(ExternalToolImporter.get_rootdir("DeLinker"))
DEEP_MODEL_FNAME = os.path.join(DELINKER_ROOT, "models/pretrained_DeLinker_model.pickle")
PAINS_DATA_FNAME = os.path.join(DELINKER_ROOT, "analysis/wehi_pains.csv")

class DeLinkerWrapper():

    DUMMY_SYMBOL = "*"
    VECTOR_SYMBOL = "*"
    def __init__(self, number_of_generation_per_valid=50, n_atomPairs_attemps=3, n_cores=1, gpu_id=None,
                 interactive=False, random_seed=None):

        self.n_atomPairs_attemps = n_atomPairs_attemps
        self.number_of_generation_per_valid = number_of_generation_per_valid
        self.n_cores= n_cores
        self.gpu_id= gpu_id
        self.interactive= interactive
        self.random_seed = -1 if not random_seed else random_seed

        (self.DeLinker_test, self.frag_utils, self.rdkit_conf_parallel,
         self.data_prepare_data, self.example_utils) = ExternalToolImporter.import_tool("DeLinker",
                                                                          ["DeLinker_test", "frag_utils",
                                                                           "rdkit_conf_parallel", "data.prepare_data",
                                                                           "example_utils"])

    def _load_andOr_prepare_mol(self, fnameOrMol, removeDummy=True, removeHs= True):
        #TODO: check if removeDummy interferes with covalent warheads
        if isinstance(fnameOrMol, str):
            mol = Chem.MolFromMolFile(fnameOrMol)
        else:
            mol = fnameOrMol
        dummies = []
        if removeDummy:
          edit_mol = Chem.RWMol(mol)
          for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetSymbol() == DeLinkerWrapper.DUMMY_SYMBOL:
              edit_mol.RemoveAtom(i)
              dummies.append( i )

            mol= edit_mol.GetMol()
        if removeHs:
          mol= Chem.RemoveHs( mol )
        mol.dummies = dummies
        return mol


    def _get_non_available_atoms(self, mol):

        for i in range( mol.GetNumAtoms() ):
            atom = mol.GetAtomWithIdx(i)
            if atom.GetImplicitValence() <1:
                yield i

    def _pick_closest_atoms(self, mol1, mol2, atom_idx1=None, atom_idx2=None, n_to_retrieve=None, cumm_idxs=True):
        '''

        :param mol1:
        :param mol2:
        :param atom_idx1:
        :param atom_idx2:
        :param n_to_retrieve:
        :param cumm_idxs: If true, indices of second molecule will start at N_ATOMS_1 position
        :return:
        '''

        #TODO: avoid calculations for valence incompatible results
        if atom_idx1 and atom_idx2:
            return atom_idx1, atom_idx2

        # for c in mol1.GetConformers():
        #     print( c.GetPositions( ))

        combo = Chem.CombineMols(mol1, mol2)
        distance_matrix = Chem.Get3DDistanceMatrix(combo)
        n_atoms1 = mol1.GetNumAtoms()
        n_atoms2 = mol2.GetNumAtoms()

        if self.interactive:
            mol_with_idxs = self.example_utils.mol_with_atom_index(combo)
            plt.imshow(Draw.MolsToGridImage([mol_with_idxs], molsPerRow=1))
            plt.show()

        distance_matrix = distance_matrix[:n_atoms1, n_atoms1:].copy()

        if atom_idx1:
            distance_matrix_tmp = np.ones_like(distance_matrix) * sys.maxsize
            distance_matrix_tmp[atom_idx1, :] = distance_matrix
            distance_matrix = distance_matrix_tmp

        if atom_idx2:
            distance_matrix_tmp = np.ones_like(distance_matrix) * sys.maxsize
            distance_matrix_tmp[:, atom_idx1] = distance_matrix
            distance_matrix = distance_matrix_tmp

        bad_atoms1 = list(self._get_non_available_atoms(mol1))
        bad_atoms2 = list(self._get_non_available_atoms(mol2))

        distance_matrix[bad_atoms1, :] =  sys.maxsize
        distance_matrix[:, bad_atoms2] =  sys.maxsize

        unravel_idxs= np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)

        selected_pairs=[]
        k = 0
        if not n_to_retrieve:
            n_to_retrieve = len(unravel_idxs[0])

        dummies_i = set(mol1.dummies)
        dummies_j = set(mol2.dummies)

        for i,j in zip(* unravel_idxs ):
            if i in dummies_i or j in dummies_j : continue
            k += 1
            if k > n_to_retrieve:
                break
            selected_pairs.append( ((int(i), int(j+ int(cumm_idxs)*n_atoms1)), distance_matrix[i,j]) )

        return selected_pairs

    def _link_molecule_pairs_given_atom_idxs(self, frag_1_mol: Chem.Mol, frag_2_mol: Chem.Mol, atom_idx_1:int, atom_idx_2:int,
                                                    min_atoms: int=2, max_atoms: int=6) -> List[Chem.Mol]:

        # RDLogger.DisableLog('rdApp.info')

        combo_no_exit = Chem.CombineMols(frag_1_mol, frag_2_mol)
        combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles(DeLinkerWrapper.VECTOR_SYMBOL+"."+DeLinkerWrapper.VECTOR_SYMBOL))  #TODO: Check if *.* can be used due to DUMMY symbol

        combo_2d = Chem.Mol(combo)
        AllChem.Compute2DCoords(combo_2d)
        # num_heavy_atoms = combo.GetNumHeavyAtoms()

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
        du = Chem.MolFromSmiles(DeLinkerWrapper.VECTOR_SYMBOL)
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
        dist, ang = self.frag_utils.compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))

        if self.interactive:
            print(Chem.MolToSmiles(mol_to_link), dist, ang)

        cwdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdirname:
            os.chdir(tmpdirname)
            data_path = "./fragments_test_data.txt"
            with open(data_path, 'w') as f:
                f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
            raw_data = self.data_prepare_data.read_file(data_path)
            self.data_prepare_data.preprocess(raw_data, "zinc", "fragments_test", True)

            if not self.gpu_id:
                os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
            else:
                os.environ['CUDA_VISIBLE_DEVICES'] = str(self.gpu_id)
                os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = "true"

            os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

            from tensorflow.python.util import deprecation
            deprecation._PRINT_DEPRECATION_WARNINGS = False

            # Arguments for DeLinker_to_remove
            args = defaultdict(None)
            args['--dataset'] = 'zinc'
            args['--config'] = '{"generation": true, \
                                               "batch_size": 1, \
                                               "number_of_generation_per_valid": %d, \
                                               "min_atoms": %d, "max_atoms": %d, \
                                               "train_file": "molecules_fragments_test.json", \
                                               "valid_file": "molecules_fragments_test.json", \
                                               "output_name": "DeLinker_example_generation.smi"}'%(self.number_of_generation_per_valid,
                                                                                                   min_atoms, max_atoms )
            args['--freeze-graph-model'] = False
            args['--restore'] = DEEP_MODEL_FNAME

            # Setup model and generate molecules

            model = self.DeLinker_test.DenseGGNNChemModel(args)
            generated_smiles = model.generate_new_graphs(model.valid_data, return_list=True)
            del model

        os.chdir(cwdir)
        in_mols = [smi[1] for smi in generated_smiles]
        frag_mols = [smi[0] for smi in generated_smiles]
        gen_mols = [smi[2] for smi in generated_smiles]

        du = Chem.MolFromSmiles(DeLinkerWrapper.VECTOR_SYMBOL)
        # for smi in frag_mols: print(  smi )
        clean_frags = [Chem.MolToSmiles(Chem.RemoveHs(
            AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi), du, Chem.MolFromSmiles('[H]'), True)[0])) for smi in
                       frag_mols]

        if self.interactive:
            print("Deep learning done!")

        # Check valid
        results = []
        for in_mol, frag_mol, gen_mol, clean_frag in zip(in_mols, frag_mols, gen_mols, clean_frags):
            if len(Chem.MolFromSmiles(gen_mol).GetSubstructMatch(Chem.MolFromSmiles(clean_frag))) > 0:
                results.append([in_mol, frag_mol, gen_mol, clean_frag])

        if self.interactive:
            print("Number of generated SMILES: \t%d" % len(generated_smiles))
            print("Number of valid SMILES: \t%d" % len(results))
            print("%% Valid: \t\t\t%.2f%%" % (len(results) / len(generated_smiles) * 100))


        linkers = list(map(self.frag_utils.get_linker, [Chem.MolFromSmiles(m[2]) for m in results],
                           [Chem.MolFromSmiles(m[3]) for m in results], [m[1] for m in results]
                           ))
        # Standardise linkers
        for i, linker in enumerate(linkers):
            if linker == "":
                continue
            linker = Chem.MolFromSmiles(re.sub('[0-9]+'+re.escape(DeLinkerWrapper.VECTOR_SYMBOL), DeLinkerWrapper.VECTOR_SYMBOL, linker))
            Chem.rdmolops.RemoveStereochemistry(linker)

            for i in range(linker.GetNumAtoms()):
                atom = linker.GetAtomWithIdx(i)
                atom.SetProp("is_linker", "True")

            linkers[i] = Chem.MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker))
        # Update results
        for i in range(len(results)):
            results[i].append(linkers[i])

        # Create dictionary of results
        results_dict = {}
        for res in results:
            if res[0] + '.' + res[
                1] in results_dict:  # Unique identifier - starting fragments and original molecule
                results_dict[res[0] + '.' + res[1]].append(tuple(res))
            else:
                results_dict[res[0] + '.' + res[1]] = [tuple(res)]

        if len(results_dict) == 0:
            return []

        # Check uniqueness
        if self.interactive:
            print("Unique molecules: %.2f%%" % (self.frag_utils.unique(results_dict.values()) * 100))

        # Check if molecules pass 2D filters
        filters_2d = self.frag_utils.calc_filters_2d_dataset(results,
                                                        pains_smarts_loc=PAINS_DATA_FNAME,
                                                        n_cores= self.n_cores)

        results_filt = []
        for res, filt in zip(results, filters_2d):
            if filt[0] and filt[1] and filt[2]:
                results_filt.append(res)

        if self.interactive:
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
            print("Number molecules passing 2D filters:\t\t%d" % len(results_filt))

        results_filt_unique = self.example_utils.unique_mols(results_filt)
        new_gen_mols_smi = [x[2] for x in results_filt_unique]

        if self.interactive:
            print("Number unique molecules passing 2D filters:\t%d" % len(results_filt_unique))

        mols_with_conf = []
        for smi in new_gen_mols_smi:
            mol = Chem.MolFromSmiles(smi)
            mol = self.place_linked_mol( mol, combo_no_exit)
            if mol:
                mols_with_conf.append(  mol )

        return mols_with_conf

    def place_linked_mol(self, mol, combined_mol_ref,):
        try:

            mol = AllChem.ConstrainedEmbed( Chem.Mol(mol), combined_mol_ref, self.random_seed)
            # print( mol.GetProp("EmbedRMS"))
            # print( AllChem.UFFGetMoleculeForceField(mol).CalcEnergy() )
            # print( AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()/mol.GetNumBonds() )
            # print( mol )
            # import matplotlib.pyplot as plt
            # from rdkit.Chem import Draw
            # plt.imshow(Draw.MolsToGridImage([mol, combined_mol_ref], molsPerRow=1)); plt.show()
            return mol
        except ValueError:
            return None


    def _estimate_num_atoms(self, distance, min_num=1, max_factor=3):
        required_atoms = round(distance/1.2)
        min_atoms = max(min_num, int(required_atoms))
        max_atoms = int(max_factor*min_atoms)
        return min_atoms, max_atoms

    def link_molecule_pair(self, fnameOrMol_1: Union[str,Chem.Mol], fnameOrMol_2:Union[str, Chem.Mol], n_attempts=None) -> List[Chem.Mol]:
        '''
        :param fnameOrMol_1:
        :param fnameOrMol_2:
        :param atom_idx_1:
        :param atom_idx_2:
        :return:
        '''
        if not n_attempts:
            n_attempts = self.n_atomPairs_attemps

        frag_1_mol = self._load_andOr_prepare_mol(fnameOrMol_1)
        frag_2_mol = self._load_andOr_prepare_mol(fnameOrMol_2)

        linked_molecules = []
        for (atom_idx_1, atom_idx_2), distance in self._pick_closest_atoms(frag_1_mol, frag_2_mol):

            min_atoms, max_atoms = self._estimate_num_atoms(distance)
            linked_molecules += self._link_molecule_pairs_given_atom_idxs( frag_1_mol, frag_2_mol, atom_idx_1, atom_idx_2,
                                                                           min_atoms, max_atoms)

            #Try with aromatic for min_atoms:
            min_atoms += 6
            max_atoms = max(min_atoms+1, max_atoms)
            linked_molecules += self._link_molecule_pairs_given_atom_idxs( frag_1_mol, frag_2_mol, atom_idx_1, atom_idx_2,
                                                                           min_atoms, max_atoms)


            n_attempts-=1
            if n_attempts<1:
                break

        linked_molecules = ( elem[0] for elem in self.example_utils.unique_mols( map(lambda mol: (mol,),linked_molecules)) )
        linked_molecules = list( filter(None.__ne__, linked_molecules ))
        return linked_molecules


def example():

    frag_binaries= [b"\xef\xbe\xad\xde\x00\x00\x00\x00\x0c\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x11\x00\x00\x00\x11\x00\x00\x00\x80\x01\x07\x00(\x00\x00\x00\x03\x03\x06\x00`\x00\x00\x00\x01\x03\x08\x00(\x00\x00\x00\x03\x02\x06\x00`\x00\x00\x00\x03\x01\x07\x00(\x00\x00\x00\x03\x03\x08\x00(\x00\x00\x00\x03\x02\x06\x00`\x00\x00\x00\x01\x03\x06\x00`\x00\x00\x00\x01\x03\x06\x00(\x00\x00\x00\x03\x04\x06\x00`\x00\x00\x00\x03\x01\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x06\x00(\x00\x00\x00\x03\x04\x06\x00`\x00\x00\x00\x02\x02\x000(\x00\x00\x00\x00\x01\x19\x08\x00\x00\x00AtomNull(\x1c\x00+\x0b\x00\x03\x00\x03\x01\x00\x03\x06\x00\x00\x07\x00\x00\x08 \x02\x08(\x02\t\x08\x00\t\n\x00\x04\x0b\x00\n\x0b\x00\x04\x0c\x00\t\r\x00\x0c\r\x00\x04\x0e \x05\x0e(\x02\x0e\x0f\x00\x0f\x10\x00\x14\x01\x06\x04\x0c\r\t\n\x0b\x17\x01\x00\x00\x00\x01\x00\x00\x00\x00\x11)\\\xdf@\xf8SS\xc01\x08\xd4A\xf8S\x0fAff\x96\xc0\x1f\x85\xd8A33\x9f@\x1f\x85K\xc0\xfc\xa9\xcbA\x7fj\xf0@;\xdf\x8b\xc0\x87\x16\xdbA\x9a\x99\xdd@\xf4\xfd\x8c\xc0\x8f\xc2\xaaA\xfa~\xb2@\xee|\x9f\xc0%\x06\x9dA1\x08\xd4@\x02+\xb3\xc0\xecQ\xdbA;\xdf\xef@5^\xfa\xbf\x0e-\xd6A\x87\x16\xc5@\xac\x1cb\xc0\x08\xac\xcbA\xb6\xf3\xd9@\xf8S\x83\xc0!\xb0\xc1A\n\xd7\xbb@\xdfO\xa1\xc0\xa6\x9b\xbbA%\x06\xd1@}?\xb1\xc0\x1f\x85\xb1A\x1f\x85\xfb@\x8bl_\xc0b\x10\xb0A'1\xe8@\xcb\xa1=\xc0\xc7K\xbaA\x85\xeb\xcd@\n\xd7\x87\xc0\x83\xc0\xa0A\xee|\xdf@\x02+G\xc0\xe9&\x9aA\x0c\x02\tA5^b\xc0\x17\xd9\x93A\x16", b'\xef\xbe\xad\xde\x00\x00\x00\x00\x0c\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x11\x00\x00\x00\x80\x01\x06\x00`\x00\x00\x00\x01\x03\x06\x00(\x00\x00\x00\x03\x04\x08\x00(\x00\x00\x00\x03\x02\x07\x00h\x00\x00\x00\x03\x02\x01\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\x07@8\x00\x00\x00\x03\x01\x03\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\x06@(\x00\x00\x00\x03\x04\t\x00 \x00\x00\x00\x01\x0b\x00\x01\x00\x01\x02(\x02\x01\x03 \x03\x04\x00\x04\x05\x00\x05\x06\x00\x06\x07h\x0c\x07\x08h\x0c\x08\th\x0c\t\nh\x0c\n\x0bh\x0c\x0b\x0ch\x0c\x0c\rh\x0c\x06\x0eh\x0c\t\x0eh\x0c\r\x0eh\x0c\x0c\x0f\x00\x14\x02\x05\x06\x07\x08\t\x0e\x06\n\x0b\x0c\r\x0e\t\x17\x01\x00\x00\x00\x01\x00\x00\x00\x00\x10b\x10\x10A\x04V\xc6@{\x14\xb5A\\\x8f$A\x9c\xc4\xbc@q=\xbbA\x1b/%A\xaa\xf1\xc2@P\x8d\xc4A\xf4\xfd6Aw\xbe\xab@\x81\x95\xb5A\x17\xd9JA\xe9&\xa1@\xdb\xf9\xbaAy\xe9PAy\xe9f@\xd9\xce\xb7A\xb8\x1eAA7\x89!@\x12\x83\xbbA\xb8\x1e/A\x91\xed,@\xe5\xd0\xc1A\x14\xae%A33\xc3?\xd1"\xc3A\xb4\xc80A\xcd\xcc\x0c?y\xe9\xbdA\xe9&-A\xfe\xd4X\xbf\xaa\xf1\xbcA\xaeG;A\xc5 \xd0\xbf\x00\x00\xb7A\xaa\xf1LAX9\x84\xbf1\x08\xb2A\xfa~PAD\x8b\xac>=\n\xb3A\xc5 BA\x85\xeb\x91?o\x12\xb9A\x17\xd9ZA\xb8\x1e\xe5\xbf33\xacA\x16']
    frag_mols = list(map(Chem.Mol, frag_binaries))
    dlw= DeLinkerWrapper(n_cores=4, number_of_generation_per_valid=10, n_atomPairs_attemps=1)
    sp= dlw.link_molecule_pair(dlw._load_andOr_prepare_mol(frag_mols[0]), dlw._load_andOr_prepare_mol(frag_mols[1]))
    print( sp )

    view_batch_size=4
    mols = []
    for i, mol in enumerate(sp):
        mols.append(mol)
        Chem.MolToMolFile(mol, os.path.expanduser("~/tmp/mols/%d.mol"%i))
        if (i+1)%view_batch_size==0:
            plt.imshow(Draw.MolsToGridImage(mols, molsPerRow=1))
            plt.show()
            mols = []

if __name__ == "__main__":
    example()

    '''

python -m fragmenstein.external.DeLinker.DeLinkerWrapper

    '''