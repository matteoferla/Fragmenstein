
import os
from typing import List

from rdkit import Chem
from rdkit import RDLogger
from timeout_decorator import timeout_decorator

from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_abstract import ErrorInComputation, CombineMerge_Base
from fragmenstein.protocols.steps.minimizePdbComplex_pyrosetta import MinimizePDBComplex_pyrosetta
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.scoring._fragmenstein_scoring import _FragmensteinScorer
from fragmenstein.utils.config_manager import ConfigManager


class CombineMerge_DeLinkerDefault( CombineMerge_Base  ):



    @staticmethod
    def get_examples_combine_params():
        data_dir = os.path.abspath(os.path.join(Xchem_info.examples_dir, "hit_mols"))
        fnames = [os.path.join(data_dir, "Mpro-x0678.mol"), os.path.join(data_dir, "Mpro-x0434.mol")]
        list_of_fragments = [Compound.MolFromFile(fname, fname.split("-")[-1].split(".")[0]) for fname in fnames]

        list_of_fragments= [ Compound.DeleteSubstructs( mol, Chem.MolFromSmiles("O=CN"))  for mol in list_of_fragments]
        list_of_fragments[0] = Compound.GetMolFrags(list_of_fragments[0])[0]
        list_of_fragments[1] = Compound.GetMolFrags(list_of_fragments[1])[1]
        list_of_fragments = [list_of_fragments]
        return dict(
            list_of_fragments= list_of_fragments
        )


    def __init__(self, random_seed=None, gpu_id=None, number_of_generation_per_valid=10, n_atomPairs_attemps=3, *args, **kwargs):

        super().__init__( *args, **kwargs)

        self.n_atomPairs_attemps = n_atomPairs_attemps
        self.number_of_generation_per_valid = number_of_generation_per_valid
        self.gpu_id = gpu_id
        self.random_seed = random_seed

        self._deLinker = None

    @property
    def deLinker(self):
        if not self._deLinker:
            from fragmenstein.external.DeLinker.DeLinkerWrapper import DeLinkerWrapper
            self._deLinker = DeLinkerWrapper( number_of_generation_per_valid=self.number_of_generation_per_valid,
                             n_atomPairs_attemps=self.n_atomPairs_attemps, n_cores=1, gpu_id=self.gpu_id,
                             interactive=False, random_seed= self.random_seed)
        return self._deLinker

    def tryOneGeneric(self, merge_id, templateFname, fragments: List[Compound], wdir, *args, **kwargs):

        assert len(fragments) == 2, "Error, DeLinker only works for pairs of compounds"

        final_outdir =  os.path.join(wdir, merge_id)
        if not os.path.exists(final_outdir):
            os.mkdir( final_outdir ) #Required, as here is where _tryOneGeneric will save checkpoint

        fragments  = list( fragments )
        for i, frag in enumerate(fragments):
            try:
                Chem.MolToMolFile( frag, os.path.join(final_outdir, "frag_%d.mol"%i))
            except Chem.rdchem.MolSanitizeException:
                pass


        RDLogger.DisableLog('rdApp.warning')
        proposed_mols = self.deLinker.link_molecule_pair(*fragments)
        RDLogger.EnableLog('rdApp.warning')

        # proposed_mols = proposed_mols[:20]
        # print( [Chem.MolToSmiles(mol) for mol in proposed_mols] )
        # import matplotlib.pyplot as plt
        # from rdkit.Chem import Draw
        # for i in range(len(proposed_mols)):
        #     plt.imshow(Draw.MolsToGridImage([proposed_mols[i]], molsPerRow=1)); plt.show()

        # proposed_mols = list(map(Chem.Mol, type(self).example_delkinker))


        placed_results = []

        w_DeLinker = Chem.SDWriter( os.path.join(final_outdir, "delinker_mols.sdf"))

        minimizer = MinimizePDBComplex_pyrosetta(templateFname, atom_constrain_filter=lambda atom: atom.HasProp("is_original_atom"))

        def minimizeMol(molId, mol ):

            try:
                mol_metadata_dict = minimizer.minimize(mol, molId=molId, outdir=final_outdir, reference_fragments=fragments)
            except (RuntimeError, TimeoutError):
                return None

            if mol_metadata_dict is None:
                return None
            else:
                mol, metadata_dict = mol_metadata_dict

            metadata_dict = _FragmensteinScorer.old_scoring_fun(metadata_dict)[-1]
            metadata_dict = _FragmensteinScorer.new_scoring_fun(metadata_dict)[-1]

            metadata_dict["fragments"] = sorted(set([ frag.primitiveId for frag in fragments]))
            metadata_dict["ref_pdb"] = templateFname

            generated_molecule = Compound( mol, molId=molId, parents= fragments)
            generated_molecule.ref_pdb = templateFname
            generated_molecule.metadata = metadata_dict
            generated_molecule.ref_molIds =  metadata_dict["fragments"]

            return  generated_molecule


        for i, proposal in enumerate(proposed_mols):
            w_DeLinker.write(proposal)
            w_DeLinker.flush()
            minimized_mol = minimizeMol( merge_id+"_"+str(i), proposal)
            if minimized_mol:
                placed_results +=  [minimized_mol ]

            # ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])
            # Victor.work_path = wdir
            # v = Victor(fragments, pdb_filename=templateFname)
            # smi = Chem.MolToSmiles( proposal )
            # v.place(smi)#, merging_mode="full")


        w_DeLinker.close()

        if len(placed_results)>0:
            w_placed = Chem.SDWriter(os.path.join(final_outdir, "deLinker_minimised_mols.sdf"))
            for mol in placed_results:
                w_placed.write(mol)

            w_placed.close()

        print(len(placed_results)) #; input("enter")
        return placed_results


def test_applyCombine():
    init_params = CombineMerge_DeLinkerDefault.get_examples_init_params()
    init_params["number_of_generation_per_valid"] = 2
    combiner = CombineMerge_DeLinkerDefault( **init_params, use_dask=True, gpu_id=0)
    results = combiner.applyCombine( **CombineMerge_DeLinkerDefault.get_examples_combine_params())
    print("RESULTS applyCombine:")
    print( results)



if __name__ == "__main__":

    print("trying combine")
    test_applyCombine()

    '''

python -m fragmenstein.protocols.steps.combineMerge_DeLinkerDefault

    '''