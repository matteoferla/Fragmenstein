import os
import pickle
import re
import tempfile
import threading
from abc import abstractmethod, ABC
from collections import defaultdict
from itertools import chain

import shutil
from zipfile import ZipFile

import dirsync
import time

import logging
from rdkit import Chem

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.protocols.dataModel.compound import Compound
from fragmenstein.protocols.steps.combineMerge_abstract import CombineMerge_Base, ErrorInComputation
from fragmenstein.protocols.steps.hitsPreprocess_base import HitsPreprocess_base
from fragmenstein.protocols.steps.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.protocols.steps.score_combinedDefault import Score_CombinedDefault
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.utils.io_utils import apply_func_to_files

RANDOM_SEED = 121

class Protocol_mergeCombineBase(ABC):

    MERGES_SUBDIR="merges"
    SCORES_SUBDIR="scoring"
    def __init__(self, data_root_dir, output_dir, merging_mode=None, filter_out_by_num_inspirational_frags=-1,
                 use_unminimized_for_ref_pdb_metadata=False, template_pattern=None, verbose=True, *args, **kwargs):

        self.data_root_dir = os.path.expanduser(data_root_dir)

        self.output_dir = os.path.expanduser(output_dir)
        self.wdir_enumeration =  os.path.join(self.output_dir, Protocol_mergeCombineBase.MERGES_SUBDIR)
        self.wdir_scoring = os.path.join(self.output_dir, Protocol_mergeCombineBase.SCORES_SUBDIR)
        self.merging_mode = merging_mode
        self.filter_out_by_num_inspirational_frags = filter_out_by_num_inspirational_frags
        self.use_unminimized_for_ref_pdb_metadata = use_unminimized_for_ref_pdb_metadata
        self.template_pattern = template_pattern
        self.verbose = verbose
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
    def hit_ids(self, values):
        raise NotImplementedError()

    @property
    def data_loader(self):
        if not self._loader:
            extra_params = LoadInput_XchemDefault.default_params_xchem()
            if self.template_pattern is not None:
                extra_params["unboundPdb_id_pattern"] = self.template_pattern
            self._loader = LoadInput_XchemDefault(self.data_root_dir, **extra_params )
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
                all_smiles, all_hit_ids = zip(*self.smile_fragIds_list) #TODO: Make this an abstract property
                all_hit_ids = sorted(set(chain.from_iterable(all_hit_ids)))
                self._ref_hits_for_scoring = all_hit_ids
        return self._ref_hits_for_scoring


    def load_avilable_results(self):

        def readFname(fname):
            with open(fname, "rb") as f:
                result = pickle.load(f)
                if isinstance(result, ErrorInComputation):
                    return [None]
                else:
                    return result

        results = apply_func_to_files(self.wdir_enumeration, ".*"+re.escape(CombineMerge_Base.RESULT_PICKLE_TEMPLATE%"")
                                      , readFname)

        results = filter(None.__ne__,
                         chain.from_iterable( results) )
        if self.max_attemps: #TODO: check if available in placeFragmenstein #TODO: Make this an abstract property
            results = HitsPreprocess_base.take_random_from_iterator(results, self.max_attemps)
        return list(results)

    def score_results(self, results):
        print("Scoring", flush=True)
        assert len(results)>0, "Error, no valid results were obtained."

        proposed_mols = []
        already_available = set([])

        for  compound in results:
            smi = Chem.MolToSmiles(compound )
            if len(compound.getFragIds())< self.filter_out_by_num_inspirational_frags: continue
            if smi in already_available: continue # Heuristic filter for uniqueness
            already_available.add(smi)
            #Get rid of bad stuff
            compound = Compound.SanitizeMol(compound)
            proposed_mols.append(compound)
        assert  proposed_mols != [], "Error, no valid molecules for scoring were found"

        scorer = Score_CombinedDefault(fragments_dir=self.data_root_dir, to_score_dir=self.wdir_enumeration,
                                       selected_fragment_ids= self.ref_hits_for_scoring,
                                       working_dir = self.wdir_scoring,
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

        compoundName_to_pdbName = []
        with ZipFile(os.path.join(self.output_dir, 'bound_minimized_pdbs.zip'), 'w') as f:
            for compound in scored_mols:
                original_name = compound.GetProp("_Name")
                compound.SetProp("original_name", original_name )
                simplified_name =  get_simplified_mol_name( compound.molId)
                compound.SetProp("_Name",simplified_name)
                minimized_pdb = os.path.join(self.wdir_enumeration, original_name, Xchem_info.predicted_boundPdb_template%original_name)
                if not os.path.isfile(minimized_pdb):
                    subdir = "_".join(original_name.split("_")[:-1]) #for delinker
                    minimized_pdb = os.path.join(self.wdir_enumeration, subdir,
                                                    Xchem_info.predicted_boundPdb_template % original_name)

                minimized_pdbBasename = Xchem_info.predicted_boundPdb_template%simplified_name
                f.write(minimized_pdb, minimized_pdbBasename)

                if self.use_unminimized_for_ref_pdb_metadata:
                    compound.SetProp("ref_pdb", get_xchem_template_name(compound.ref_pdb))
                else:
                    compound.SetProp("ref_pdb", minimized_pdbBasename)

                compoundName_to_pdbName.append( "%s\t%s\t%s"%( simplified_name, minimized_pdbBasename, original_name))
            f.writestr("info/compoundName_to_pdbName.tab", "\n".join(compoundName_to_pdbName))

        if hasattr(self, "template_xchemId"):
            template_xchemId = self.template_xchemId
        else:
            template_xchemId = None

        frag_writer = FragalysisFormater(ref_pdb_xchemId=template_xchemId, addtitional_fields=['atomic_mapping']) #TODO: store original name and original template
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

        in_dir = os.path.abspath(os.path.expanduser(kwargs["data_root_dir"]))
        out_dir = os.path.abspath(os.path.expanduser(kwargs["output_dir"]))

        scores = None

        existing_outdir = True
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            existing_outdir = False
        if "working_dir" in kwargs and  kwargs["working_dir"] is not None:
            wdir = kwargs["working_dir"]
        else:
            wdir = tempfile.gettempdir()
        only_enumeration = ("skip_scoring" in kwargs and kwargs["skip_scoring"])
        only_evaluation = ("skip_enumeration_and_score_available" in kwargs and  kwargs["skip_enumeration_and_score_available"])

        assert (only_enumeration and only_evaluation) == False, "Error, skip_scoring and skip_enumeration_and_score_available are incompatible arguments"

        keep_sync = True
        def syncronizer(new_out_dir, time_sleep):
            while keep_sync:
                time.sleep(time_sleep)
                dirsync.sync(new_out_dir, out_dir, 'sync', verbose=False, logger=logging.getLogger('dummy'))

        with tempfile.TemporaryDirectory(dir="/dev/shm") as tmp_indir, \
             tempfile.TemporaryDirectory(dir=wdir) as tmp_outdir: #TODO: tmp_outdir should not be temporary if working_dir provided
            new_in_dir = os.path.join(tmp_indir, os.path.basename(in_dir))
            # print("Copying data to working directories (%s, %s)..."%(new_in_dir, tmp_outdir), end=" ", flush=True)
            if not os.path.exists(new_in_dir):
                os.mkdir(new_in_dir)
            dirsync.sync(in_dir, new_in_dir, 'sync', verbose=False, logger=logging.getLogger('dummy'))
            # shutil.copytree(in_dir, new_in_dir)
            kwargs["data_root_dir"] = new_in_dir


            new_out_dir = os.path.join(tmp_outdir, os.path.basename( os.path.abspath(out_dir)))
            os.mkdir(new_out_dir)
            if existing_outdir: #create links to files that could be readed
                new_scoring_dir = os.path.join(new_out_dir, cls.SCORES_SUBDIR)
                os.mkdir(new_scoring_dir)
                old_scoring_dir = os.path.join(out_dir, cls.SCORES_SUBDIR)
                if os.path.exists(old_scoring_dir):
                    for name in os.listdir(old_scoring_dir):
                        os.symlink(os.path.join(old_scoring_dir, name), os.path.join(new_scoring_dir, name))
                else:
                    os.mkdir( old_scoring_dir )

                new_merges_dir = os.path.join(new_out_dir, cls.MERGES_SUBDIR)
                old_merges_dir = os.path.join(out_dir, cls.MERGES_SUBDIR)
                if os.path.exists(old_merges_dir):
                    if  only_evaluation:
                        os.symlink(os.path.join(out_dir, cls.MERGES_SUBDIR), new_merges_dir)
                    else:
                        os.mkdir(new_merges_dir)
                        for name in os.listdir(old_merges_dir):
                            os.symlink(os.path.join(old_merges_dir, name), os.path.join(new_merges_dir, name))
                else:
                    os.mkdir(new_merges_dir)

            kwargs["output_dir"] = new_out_dir
            syncronizerThr = threading.Thread(target=syncronizer, args=(new_out_dir, 30))
            syncronizerThr.start()

            protocol = cls(*args, **kwargs)
            protocol.initialize(*args, **kwargs)

            if only_evaluation:
                results = protocol.load_avilable_results()
            else:
                results = protocol.compute()

            if not only_enumeration:
                scores = protocol.score_results(results)

            keep_sync = False
            syncronizerThr.join()
            dirsync.sync(new_out_dir, out_dir, 'sync', verbose=False )

        return scores


    @classmethod
    def generateCmdParser(cls, prog, description):
        from fragmenstein.utils.cmd_parser import ArgumentParser
        parser = ArgumentParser(prog=prog, description=description)


        parser.add_argument("-o", "--output_dir", type=str, help="The directory where results will be saved",
                            required=True)

        parser.add_argument("-d", "--templates_dir", type=str,
                            help="The directory where templates would be look for if not --template or", required=False)

        parser.add_argument("--template_pattern", type=str,
                            help="The regex pattern of a template to find them in --templates_dir. "
                                 "Default: '%(default)s'", required=False, default= Xchem_info.unboundPdb_id_pattern)

        parser.add_argument( "--skip_scoring", action="store_true",
                            help="Do not score enumerated molecules")

        parser.add_argument( "--skip_enumeration_and_score_available", action="store_true",
                            help="Do not propose more molecules and only score them ")

        parser.add_argument( "--filter_out_by_num_inspirational_frags", type=int, default=-1,
                            help="Discard compounds that only match to NUM number of one inspirational hit. Set to -1 to ignore"
                                 "this filter. Default %(default)s")

        parser.add_argument( "--working_dir", type=str, default=None, help="Directory where results are computed before"
                                                                           " being synchronized to output_dir")

        parser.add_argument( "--use_unminimized_for_ref_pdb_metadata",action="store_true",
                            help="Use the pre-minimized pdb fname for the ref_pdb medtadata field. Default %(default)s")
        return parser