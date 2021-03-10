import csv
import json
import logging
import os
import itertools

from abc import ABC, abstractmethod, abstractproperty
import dask.bag as DB

from dask.distributed import progress, as_completed

from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from rdkit import Chem

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.scoring.scorer_labels import SCORE_NAME_TEMPLATE, FRAGMENTS_ID, MOL_NAME_ID, checkIfNameIsScore
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.io_utils import load_mol_if_str
from fragmenstein.utils.parallel_utils import get_parallel_client

journal = logging.getLogger('Scorer')
journal.setLevel(logging.DEBUG)



computeters = []
def joblibMapFunction(cls, c_args, c_kwargs, inner_args):
    print(cls, c_args, c_kwargs, inner_args) #; input("enter")
    if len(computeters) == 0:
        computeters.append(cls(*c_args, **c_kwargs))
    computer = computeters[0]
    mol_id, (mol, frag_ids) = inner_args
    return computer.processOneMolecule(mol_id, mol, frag_ids)


class _ScorerBase(ABC):
    '''
    In order to use this _ScorerUtils_mixin, create a child class and implement the method computeScoreOneMolecule(self, mol_id, mol, frags_dict, *args, *kwargs)
    '''

    MOL_NAME_ID = MOL_NAME_ID
    FRAGMENTS_ID = FRAGMENTS_ID
    SCORE_NAME_TEMPLATE = SCORE_NAME_TEMPLATE

    @classmethod
    def parseCmd(cls, description, additional_args: List[Tuple[str,str,Dict]] = {}):
        '''

        :param description: Program description
        :param additional_args:
                                    e.g.

                                    [

                                        ("-i", "--input", {
                                          "type": str,
                                          "nargs": None,
                                          "required": True,
                                          "help": "help msg n 1"}),

                                      ("-o", "--output", {
                                          "type": str,
                                          "nargs": None,
                                          "required": True,
                                          "help": "help msg 2"}),

                                        ("-p", "--processingType", {
                                          "choices": ['wideTarget', 'tightTarget', 'highRes'],
                                          "default": 'tightTarget',
                                          "help": "help msg 3"}),
                                    ]
        :return: a dictionary with the parsed arguments
        '''

        #input_file example
        '''
        mol_id, mol_fname, fragment_ids
        x10322_0A,%(input_dir)s/Mpro-x10322_0A/Mpro-x10322_0A.mol,"x1234,x1324"
        x11810_0A,%(input_dir)s/Mpro-x11810_0A/Mpro-x11810_0A.mol,nan
        x0464_0A,%(input_dir)s/Mpro-x0464_0A/Mpro-x0464_0A.mol,x1234
        x12177_0A,%(input_dir)s/Mpro-x12177_0A/Mpro-x12177_0A.mol,nan
        '''

        import argparse
        parser = argparse.ArgumentParser(description=description)

        for args in additional_args:
            parser.add_argument(*args[:-1], **args[-1])

        def add_argument(*args, **kwargs):
            try:
                parser.add_argument(*args, **kwargs)
            except argparse.ArgumentError:
                pass

        #TODO: input should be  either a csv file or a sdf
        add_argument('-i', '--input_file', required=True, type=str, help='File containing the name and filename of the molecules to analize.'
                                                                       ' Format: csv file with the following columns: "mol_id, mol_fname, fragment_ids".'
                                                                       'If nan in fragment_ids columns, all fragments will be used for evaluation')

        add_argument('-d', '--input_dir', required=False, type=str, help='Replace "%%(input_dir)s" to INPUT_DIR for the fnames stored in mol_fnames column')

        add_argument('-f', '--fragments_dir', required=True, type=str, help='Directory with a file for each fragment to ' #This  help may be overwritten by subclasses
                                                                     'compare. Fragments can also be located in subdirectories)')

        add_argument('-p', '--fragment_id_pattern', required=False, help='Regex pattern for the fragment files.', default=r"(.+)\.mol$")

        add_argument('-o', '--table_output', required=False, default=None, type=str, help='Fname for results table')

        add_argument('-s', '--sdf_output', required=False, default=None, type=str, help='Fname for an sdf file of the evaluated molecules')

        add_argument('-w', '--working_dir', required=True, type=str, help='Directory where partial results will be saved')


        args = parser.parse_args()
        return vars( args )

    @classmethod
    def evalPipeline(cls, initiaze_parallel_execution=True):

        args = cls.parseCmd()
        journal.warning( args )

        if initiaze_parallel_execution:
            dask_client = get_parallel_client()
            print(dask_client)

        input_table = pd.read_csv( args["input_file"] )

        if "input_dir" in args:
            input_dir = {"input_dir":args["input_dir"]}
            input_table["mol_fname"]= input_table["mol_fname"].map( lambda x: x%input_dir  , na_action="ignore" )

        def processFragList(fragLists):
            if isinstance(fragLists, str):
                return fragLists.split(",")
            elif np.isnan(fragLists):
                return None
            else:
                raise ValueError("Error in file format %s"%(args["input_file"]))

        molId_to_molAndfragIds = { molId:(fname, processFragList(fragLists)) for molId,fname,fragLists in input_table.values }

        assert  args["table_output"] is not None or args["sdf_output"] is not None, "Error, at least one of the following outputs should be provided: 'table_output', 'sdf_output'"
        results_sdf_fname = args["sdf_output"]
        results_table_fname = args["table_output"]
        results = cls.computeScoreForMolecules(molId_to_molAndfragIds, results_sdf_fname= results_sdf_fname, results_table_fname= results_table_fname, **args)
        if initiaze_parallel_execution:
            dask_client.close()
        return results

    @classmethod
    def computeScoreForMolecules(cls, molId_to_molAndfragIds: Dict[str, Tuple[Chem.Mol, Union[List[str], None]]],
                                 results_sdf_fname: str = None, ref_pdb_xchemId: str=None, results_table_fname: str= None,
                                 *args, **kwargs):
        '''
        :param molId_to_molAndfragIds: dict of molecules to evaluate. mol_id -> (mol, [frag_ids]). If frag_ids is None, use all fragments
                                mol can be either a Chem.Mol or a filename
        :param results_sdf_fname: an sdf file where the molecules would be stored and the scores, fragments found and additional
                                information would be stored as molecule properties. First molecule is a dummy molecule that describe
                                the fields.

        :param ref_pdb_xchemId: string : The xchem Id for the atomic model used during computations. Usesd only to write sdf file. If None,
                                         the first fragment from each element in  molId_to_molAndfragIds would be used instead

        :param results_table_fname: string : fname where summary table will be saved.


        :return:
        '''
        computer = cls(*args, **kwargs)
        alreadyComputed_or_None = list(map(computer.loadPreviousResult, molId_to_molAndfragIds.keys()))
        not_computed_mols =  (mol_and_info for elem, mol_and_info in zip(alreadyComputed_or_None, molId_to_molAndfragIds.items()) if elem is None)

        # alreadyComputed_or_NoneD = { x["name"]: x for x in alreadyComputed_or_None}
        # input( alreadyComputed_or_NoneD["x0438-0B-x0283-0B-x0183-0B-0b0-x0183-0B-0b1-x0183-0B-0b3"] )


        def mapFunction(args):
            mol_id, (mol, frag_ids) = args
            return computer.processOneMolecule(mol_id, mol, frag_ids)

        results_new = DB.from_sequence(not_computed_mols).map(mapFunction)  # .filter(keep_fun)
        prev_results= DB.from_sequence(alreadyComputed_or_None).filter(None.__ne__)

        results_computed = DB.concat([results_new, prev_results])

        scores_ids = None
        record = None

        results_computed = list(results_computed) #We want to store it in memory as a list to return it

        for record in results_computed:
            if record is None: continue
            scores_ids = [ elem for elem in record.keys() if checkIfNameIsScore(elem) ]
            break
        assert  scores_ids is not None, "Error, not even a single molecule was scored"


        if results_sdf_fname or results_table_fname:
            md_dicts_list, mols_list = [], []
            if ref_pdb_xchemId:
                add_ref_pdb = lambda mol, frag_ids: mol
            else:
                def add_ref_pdb(mol, frag_ids):
                    mol.SetProp(FragalysisFormater.REF_PDB_FIELD , frag_ids[0].split("_")[0] )
                    return mol

            panda_rows = []
            for record in itertools.chain.from_iterable([[record], results_computed]):
                input(record)
                if record is None: continue
                panda_rows.append( [record[_ScorerBase.MOL_NAME_ID]] + [ record[score_id] for score_id in scores_ids] + [  ",".join(record[_ScorerBase.FRAGMENTS_ID]) ] )
                mol_name = record[_ScorerBase.MOL_NAME_ID]
                mol, fragments = molId_to_molAndfragIds[mol_name]
                if fragments is None:
                    fragments = record[_ScorerBase.FRAGMENTS_ID]
                if len(fragments) == 0: continue
                mol = load_mol_if_str(mol)
                if mol is None: continue

                mol = add_ref_pdb(mol, fragments)
                md_dicts_list.append( record )
                mols_list.append( mol )

            assert  len(md_dicts_list) >0, "Error, no molecules scored"
            if results_sdf_fname:
                FragalysisFormater(ref_pdb_xchemId= ref_pdb_xchemId).write_molsList_to_sdf( results_sdf_fname, mols_list, md_dicts_list)

            if results_table_fname:
                df = pd.DataFrame(panda_rows, columns=[_ScorerBase.MOL_NAME_ID]+scores_ids+ [_ScorerBase.FRAGMENTS_ID])
                df.sort_values(by=scores_ids[0], inplace=True)
                df.to_csv(results_table_fname, index=False,quoting=csv.QUOTE_NONNUMERIC)

        results_computed = { record[cls.MOL_NAME_ID]: record for record in results_computed }
        return results_computed

    def __init__(self, working_dir, verbose=False, *args, **kwargs):
        if not os.path.exists(working_dir):
            raise ValueError(("Error, working directory (%s) does not exists. If this is the first time you execute the progam, "+
                             "please, create the directory")%working_dir)
        self.working_dir = working_dir
        self.verbose = verbose
        # assert  hasattr(self, " fragments_id"), "Error,  fragments_id is a required attribute, but was not used"

    @abstractproperty
    def fragments_id(self):
        raise NotImplementedError("Error, this is base clase")

    def getCheckpointName(self, mol_or_id, wdir=None):
        if isinstance(mol_or_id, str):
            name = mol_or_id
        else:
            journal.warning("Not recommended option (hash mol) for getCheckpointName")
            name = str(abs(hash(mol_or_id)))  # TODO: This could be not unique enough
        wdir = self.working_dir if wdir is None else wdir
        return os.path.join(wdir, name + "_%s.json"%str(type(self).__name__))

    def loadPreviousResult(self,  mol_or_id ):
        checkpoint_fname = self.getCheckpointName( mol_or_id)
        try:
            with open(checkpoint_fname) as f:
                data = json.load(f)
            return data
        except FileNotFoundError:
            return None


    def processOneMolecule(self, mol_id, mol, frag_ids, *args, **kwargs):

        # print( "processing %s"%mol_id)
        result =  self.loadPreviousResult(mol_id)
        if result is not None:
            return result

        mol = load_mol_if_str(mol)
        if mol is None:
            result = None
        else:

            if frag_ids is None:
                frag_ids =  self.fragments_id

            result = self.computeScoreOneMolecule(mol_id, mol, frag_ids, *args, **kwargs)

        with open( self.getCheckpointName( mol_id), "w") as f:
            json.dump(result, f)
        return result

    @abstractmethod
    def computeScoreOneMolecule(self, mol_id, mol, frag_ids, *args, **kwargs) -> Dict[str, Tuple[float, List[str]]]:
        '''

        :param mol_id:
        :param mol: the molecule to be scored
        :param frag_ids: the xchem ids for the inspirationl hits
        :param args:
        :param kwargs:
        :return: A dict with the following structure {_ScorerBase.MOL_NAME_ID: mol_id, "score_xcos": score, _ScorerBase.FRAGMENTS_ID: list(fragment_ids) }
        '''
        raise NotImplementedError("This is an abstact method that have to be overwritten")
