import os
import re
from abc import abstractmethod, ABC
from collections import defaultdict
from itertools import chain

from rdkit import Chem

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.protocols.steps.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.protocols.steps.score_combinedDefault import Score_CombinedDefault
from fragmenstein.protocols.xchem_info import Xchem_info

RANDOM_SEED = 121

class Protocol_mergeCombineBase(ABC):

    def __init__(self, data_root_dir, output_dir, merging_mode=None, *args, **kwargs):

        self.data_root_dir = os.path.expanduser(data_root_dir)

        self.output_dir = os.path.expanduser(output_dir)
        self.wdir_fragmenstein =  os.path.join(self.output_dir, "merges")
        self.merging_mode = merging_mode

        self._loader = None
        self._fragments = None
        self._fragments_dict = None
        self._ref_hits_for_scoring = None
        self.mol_num = 0
        self._cur_molNum_scored = 0


    @property
    def cur_molNum_scored(self):
        '''
        Used for scoring only
        :return:
        '''
        val = self._cur_molNum_scored
        self._cur_molNum_scored += 1
        return val

    @property
    @abstractmethod
    def sdf_outname(self):
        raise NotImplementedError()

    @property
    @abstractmethod
    def hit_ids(self):
        raise NotImplementedError()

    @property
    def data_loader(self):
        if not self._loader:
            self._loader = LoadInput_XchemDefault(self.data_root_dir, **LoadInput_XchemDefault.default_params_xchem() )
        return self._loader

    @property
    def fragments(self):

        if not self._fragments:
            if self.hit_ids:
                self._fragments = list(filter(lambda frag: frag.molId in self.hit_ids, self.data_loader.fragments))
            else:
                self._fragments = self.data_loader.fragments
        return self._fragments

    @property
    def fragments_dict(self):

        if not self._fragments_dict:
            if self.hit_ids:
                self._fragments_dict = {fragId:frag for fragId, frag in self.data_loader.fragments_dict.items() if fragId in self.hit_ids}
            else:
                self._fragments_dict = self.data_loader.fragments_dict
        return self._fragments_dict


    @property
    def ref_hits_for_scoring(self):

        if not self._ref_hits_for_scoring:
            if hasattr(self, "hit_ids"):
                self._ref_hits_for_scoring = self.hit_ids
            else:
                all_smiles, all_hit_ids = zip(*self.smile_fragIds_list)
                all_hit_ids = sorted(set(chain.from_iterable(all_hit_ids)))
                self._ref_hits_for_scoring = all_hit_ids
        return self._ref_hits_for_scoring


    def score_results(self, results):
        print("Scoring")
        assert len(results)>0, "Error, no valid results were obtained."

        proposed_mols = []
        already_available = set([])
        for  mol in results:
            # Heuristic filter for uniqueness
            smi = Chem.MolToSmiles(mol )
            if smi in already_available: continue
            already_available.add(smi)
            proposed_mols.append(mol)

        scorer = Score_CombinedDefault(fragments_dir=self.data_root_dir, to_score_dir=self.wdir_fragmenstein,
                                       selected_fragment_ids= self.ref_hits_for_scoring,
                                       working_dir = os.path.join(self.output_dir, "scoring"),
                                       **Score_CombinedDefault.default_params_xchem())

        scored_mols = scorer.compute_scores(proposed_mols)
        print(scored_mols)

        def get_simplified_mol_name(mol_id):  # TODO: Move it within fragalysis??
            '''
            Reduce the mol_name for merges over fragments of fragments
            :param mol_id:
            :return:
            '''
            molId_additionalInfo = mol_id.split("_")
            mol_id = molId_additionalInfo[0]
            fragments = [ elem.rstrip("-") for elem in filter(lambda x: x!="", mol_id.split("x"))]
            fragIds_dict = defaultdict(lambda : defaultdict(set))
            for fragId_chainId_bitId in fragments:
                fragId_chainId_bitId = fragId_chainId_bitId.split("-")
                if len(fragId_chainId_bitId)==2:
                    fragId, chainId = fragId_chainId_bitId
                    bitId = ""
                else:
                    fragId, chainId, bitId = fragId_chainId_bitId[:3]
                fragId = "x"+ fragId
                fragIds_dict[fragId][chainId].add( bitId )

            mol_id = "%d_"%self.cur_molNum_scored
            for fragId, chain_to_bit_dict in fragIds_dict.items():
                mol_id += fragId+"_"
                for chainId, bitIds in chain_to_bit_dict.items():
                    mol_id += chainId
                    if len(bitIds)>0:
                        mol_id += "_"+ "_".join( sorted(bitIds))

                mol_id += "-"

            mol_id = mol_id.rstrip("-").rstrip("_")
            return mol_id


        def get_xchem_template_name(ref_pdb):
            ref_pdb = os.path.basename(ref_pdb)
            ref_pdb = re.match(Xchem_info.fragment_no_chain_pattern, ref_pdb).group(1)
            return ref_pdb

        for mol in scored_mols:
            mol.SetProp("original_name", mol.GetProp("_Name") )
            mol.SetProp("_Name", get_simplified_mol_name( mol.molId))
            mol.SetProp("ref_pdb", get_xchem_template_name(mol.ref_pdb))

        if hasattr(self, "template_xchemId"):
            template_xchemId = self.template_xchemId
        else:
            template_xchemId = None

        frag_writer = FragalysisFormater(ref_pdb_xchemId=template_xchemId) #TODO: store original name and original template
        frag_writer.write_molsList_to_sdf(self.sdf_outname, scored_mols) #, metadata_list)
        print("Results written to %s" % self.sdf_outname)

        return scored_mols

    @abstractmethod
    def compute(self):
        raise  NotImplementedError()

    @abstractmethod
    def initialize(self, *args, **kwargs):
        raise  NotImplementedError()

    @classmethod
    def main(cls, *args, **kwargs):

        protocol = cls(*args, **kwargs)
        protocol.initialize(*args, **kwargs)
        results = protocol.compute()

        scores = protocol.score_results(results)
        return scores
