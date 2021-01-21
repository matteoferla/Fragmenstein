import csv
import os
import warnings
import json
import dask
import numpy as np
import pandas as pd

from itertools import chain

from typing import Dict, List, Tuple, Union
from rdkit import Chem

from fragmenstein.scoring._scorer_utils import _ScorerUtils, prepare_paralell_execution


class _ScorerBase(_ScorerUtils):
    '''
    In order to use this _ScorerUtils_mixin, create a child class and implement the method computeScoreOneMolecule(self, mol_id, mol, frags_dict, *args, *kwargs)
    '''

    @classmethod
    def processOneMolecule(cls, mol_id, mol, frag_dict, *args, **kwargs):

        computer = cls(*args,**kwargs)
        result =  computer.loadPreviousResult(mol_id)
        if result is not None:
            return result

        if isinstance(mol, str):
            try:
                mol = cls.load_molecule( mol)
            except OSError:
                print("OSError", mol_id, mol, )
                mol = None
        if mol is None:
            result = None
        else:
            result = computer.computeScoreOneMolecule(mol_id, mol, frag_dict, *args, **kwargs)

        with open( computer.getCheckpointName( mol_id), "w") as f:
            json.dump(result, f)

        return result

    @classmethod
    def computeScoreForMolecules(cls, molId_to_mol_fragIds: Dict[str, Tuple[Union[Chem.Mol, str], List[str]]], frag_dict: Dict[str, Chem.Mol], results_table_fname: str= None
                                 , *args, **kwargs):
        '''
        :param molId_to_mol_fragIds: dict of molecules to evaluate. mol_id -> (mol, [frag_ids]). If frag_ids is None, use all fragments
                                mol can be either a Chem.Mol or a filename
        :param frag_dict: a dict of frag_id -> Chem.Mol to compare with bit
        :param results_table_fname: string : fname where summary table will be saved.
        :return:
        '''

        computer = cls(*args,**kwargs)


        alreadyComputed_or_None = list(map(computer.loadPreviousResult, molId_to_mol_fragIds.keys()))
        not_computed_mols =  (mol_and_info for elem, mol_and_info in zip(alreadyComputed_or_None, molId_to_mol_fragIds.items()) if elem is None)

        results_computed = []
        for mol_id, (mol, fragIds) in not_computed_mols:
            if fragIds is None:
                current_frag_dict = frag_dict
            else:
                current_frag_dict = { key: frag_dict[key] for key in fragIds}
            result = dask.delayed(cls.processOneMolecule)(mol_id, mol, current_frag_dict, *args, **kwargs)
            results_computed.append( result )
        results_computed = dask.compute(results_computed)[0]

        results_computed = chain.from_iterable( [results_computed, filter(None.__ne__, alreadyComputed_or_None) ] )


        if results_table_fname:
            panda_rows = []
            for record in results_computed:
                if record is None: continue
                panda_rows.append( [record["mol_name"], record["score"], ",".join(record["fragments"]) ] )
            df = pd.DataFrame(panda_rows, columns=["mol_name", "score", "fragments"])
            df.sort_values(by="score", inplace=True)
            df.to_csv(results_table_fname, index=False,quoting=csv.QUOTE_NONNUMERIC)

        return results_computed


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
        x10322_0A, % (input_dir)
        s / Mpro - x10322_0A / Mpro - x10322_0A.mol, x1234,x1324
        x11810_0A, % (input_dir)
        s / Mpro - x11810_0A / Mpro - x11810_0A.mol, nan
        x0464_0A, % (input_dir)
        s / Mpro - x0464_0A / Mpro - x0464_0A.mol, x1234
        x12177_0A, % (input_dir)
        s / Mpro - x12177_0A / Mpro - x12177_0A.mol, nan
        '''

        import argparse
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument('-i', '--input_file', required=True, type=str, help='File containing the name and filename of the molecules to analize.'
                                                                       ' Format: csv file with the following columns: "mol_id, mol_fname, fragment_ids".'
                                                                       'If nan in fragment_ids columns, all fragments will be used for evaluation')

        parser.add_argument('-d', '--input_dir', required=False, type=str, help='Replace "%%(input_dir)s" to INPUT_DIR for the fnames stored in mol_fnames column')

        parser.add_argument('-f', '--fragments', required=True, type=str, help='Directory with a mol file for each fragment to '
                                                                     'compare. Fragments can also be located in subdirectories)')

        parser.add_argument('-p', '--fragment_id_pattern', required=False, help='Regex pattern for the fragment files.', default=r"(.+)\.mol$")

        parser.add_argument('-o', '--output', required=True, type=str, help='Fname for results table')

        parser.add_argument('-w', '--working_dir', required=True, type=str, help='Directory where partial results will be saved')

        for kwarg in additional_args:
            parser.add_argument(*kwarg[:-1], **kwarg[-1])

        args = parser.parse_args()
        return vars( args )

    @classmethod
    def evalPipeline(cls, initiaze_parallel_execution=True):

        args = cls.parseCmd()
        print(args)

        if initiaze_parallel_execution:
            dask_client = prepare_paralell_execution()
            print(dask_client)


        fragments = cls.load_files_as_mols( args["fragments"], file_pattern=args["fragment_id_pattern"], delay_loading=False)
        fragments_dict = dict(fragments)


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
        molId_to_mol_fragIds = { molId:(fname, processFragList(fragLists)) for molId,fname,fragLists in input_table.values }

        cls.computeScoreForMolecules( molId_to_mol_fragIds, fragments_dict, results_table_fname= args["output"], wdir = args["working_dir"])

    def __init__(self, wdir):
        self.wdir = wdir

    def getCheckpointName(self, mol_or_id, wdir=None):
        if isinstance(mol_or_id, str):
            name = mol_or_id
        else:
            warnings.warn("Not recommended option (hash mol) for getCheckpointName")
            name = str(abs(hash(mol_or_id)))  # TODO: This could be not unique enough
        wdir = self.wdir if wdir is None else wdir
        return os.path.join(wdir, name + "_xcos.json")

    def loadPreviousResult(self,  mol_or_id ):
        checkpoint_fname = self.getCheckpointName( mol_or_id)
        try:
            with open(checkpoint_fname) as f:
                data = json.load(f)
            return data
        except FileNotFoundError:
            return None

    def computeScoreOneMolecule(self):
        raise NotImplementedError("This is an abstact method that have to be overwritten")
