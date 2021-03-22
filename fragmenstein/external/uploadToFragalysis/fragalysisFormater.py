import json
import logging
import os
import re
import pandas as pd

from typing import Dict, Union
from typing import List
from datetime import datetime

from collections import OrderedDict

from itertools import chain

from functools import reduce
from rdkit import Chem

from fragmenstein.scoring.scorer_labels import  checkIfNameIsScore
from fragmenstein.utils.io_utils import load_mol_if_str, apply_func_to_files, load_files_as_mols

journal = logging.getLogger('Fragalysis Formatter')
journal.setLevel(logging.DEBUG)

class FragalysisFormater():

    METADATA_FIELDS_DEFAULT_FNAME= os.path.abspath( os.path.join(__file__, os.path.pardir, "metadata_fields_example.csv"))
    REF_PDB_FIELD = "ref_pdb"
    REQUIRED_FIELDS = ["_Name", "original SMILES", REF_PDB_FIELD, "ref_mols"]
    OPTIONAL_FIELDS = ["original_name"]
    FRAGALYSIS_UPLOAD_VERSION = "ver_1.2"

    @classmethod
    def _check_property(cls, mol, prop):
        if mol.HasProp(prop) and not mol.GetProp(prop) in [None, ""]:
            return True
        else:
            return False

    @classmethod
    def _check_molecule_metadata(cls, mol, scores_names):
        for required in list(FragalysisFormater.REQUIRED_FIELDS)+list(scores_names):
            if not cls._check_property(mol, required):
                journal.error("Mol (%s) does not contain property %s\n\n"%(str(mol),required))
                return False
        return True

    @classmethod
    def _add_prop_if_not_available(cls, mol, name, val): #TODO: move this things to utils
        '''
        Adds property to a molecule if not already existing
        :param mol:
        :param name:
        :param val:
        :return:
        '''
        if not cls._check_property(mol, name):
            mol.SetProp(name, val)

    @classmethod
    def trimStr(cls, val):
        if len(val) > 50:
            new_val = val[:50]
            journal.warning(
                "Name %s is to long to be displayed in Fragalysis, it would be trimemed to %s" % (val, new_val))
            val = new_val
        return val

    @classmethod
    def trimCharFieldInMol(cls, mol, propertyName):
        val = mol.GetProp(propertyName)
        val = cls.trimStr(val)
        mol.SetProp(propertyName, val)
        return mol

    @classmethod
    def add_metadata_to_mol(cls, mol, metadata_dict, round_digits=6):
        for prop, val in metadata_dict.items():
            if checkIfNameIsScore(prop):
                val = str(round(val, ndigits=round_digits))
            elif prop == "fragments":
                prop = "ref_mols"
                val = ",".join(val)
            elif prop in  ["name", "mol_name"]:
                prop= "_Name"
                val = cls.trimStr( val )

            cls._add_prop_if_not_available(mol, prop, str(val))
        return mol

    @classmethod
    def merge_sdfs(cls, outname, *fnames):
        '''
        It does not check if the sdf files are compatible.
        :param outname:
        :param *fnames:
        :return:
        '''
        header_done = False
        n_lines_detect_header = 0
        with open(outname, "w") as f_out:
            for fname in fnames:
                with open(fname) as f:
                    header_found = False
                    for line in f:
                        if not header_found:
                            if line.startswith("$$$$"):
                                header_done=True
                                header_found=True
                            if  not header_done:
                                f_out.write( line )
                        else:
                            n_lines_detect_header += 1
                            f_out.write(line)
        assert  n_lines_detect_header>0, "Error, header was not found."

    def __init__(self, ref_pdb_xchemId=None, metadata_header=None, drop_unknown_fields=True, addtitional_fields=[]): #TOOD: Ref pdb could be many different. One per compound even

        #TODO: add option ref_pdbs_folder
        self.drop_unknown_fields = drop_unknown_fields

        if metadata_header is None:
            self.metadata_header =  self.parse_metadata_config(FragalysisFormater.METADATA_FIELDS_DEFAULT_FNAME)

        else:
            self.metadata_header = self.parse_metadata_config(metadata_header)


        self.ref_pdb_xchemId = ref_pdb_xchemId

        self.optional_fileds = list(FragalysisFormater.OPTIONAL_FIELDS)+addtitional_fields

    def parse_metadata_config(self, fname_or_iterable):

        if isinstance(fname_or_iterable, str):
            md_fields = pd.read_csv(fname_or_iterable, comment="#")
            md_fields =  OrderedDict( elem for idx, elem in md_fields.iterrows() )
        elif isinstance(fname_or_iterable, (OrderedDict, list)):
            md_fields = OrderedDict(fname_or_iterable)
        else:
            raise ValueError("fname_or_iterable must be the name of a csv file or an ordered dict ")
        md_fields["generation_date"] = datetime.today().strftime('%Y-%m-%d')
        return md_fields

    def _get_header_mol(self, missing_properties=[]):
        mol = Chem.MolFromSmiles("C")
        mol.SetProp("_Name", FragalysisFormater.FRAGALYSIS_UPLOAD_VERSION)
        for key, value in chain.from_iterable( [self.metadata_header.items(), missing_properties] ):
            mol.SetProp(key, value)
        return mol

    def _autodetect_missing_scores_in_header(self, mol):
        missing_properties = []
        for propName in mol.GetPropNames():
            if checkIfNameIsScore(propName) and not propName in self.metadata_header:
                missing_properties.append( (propName, "Score %s"%propName) )
        return missing_properties

    def write_molsList_to_sdf(self, results_sdf_fname, mol_list, metada_list: List[Dict[str, Union[str, float]]]=None):
        '''

        :param results_sdf_fname: The file where formated sdf file will be saved
        :param mol_list: a list of molecules to write. They must contain the following properties (set with mol.SetProp("name", value)
                         Properties could be, instead, contained in the parameter metadata_list, a list of dicts.

                        Mandatory:

                        _Name: The name of the molecule
                        ref_mols: Comma separated list of xchem ids for the compounds used as inspiration

                        Optional
                        ref_pdb: The xchem id available in Fragalysis for the atomic model that was used. If no provided as a parameter
                                 the default one would be used instead (self.ref_pdb_xchemId)
                        original SMILES: The smiles of the compound. If not provided, it would be automatically computed

        :return:
        '''

        if metada_list:
            assert len(metada_list) == len(mol_list), "Error, when metadata is provided as a list of dicts, the number of dicts should be" \
                                                      " equal to the number of molecules"
            for mol, md_dict in zip(mol_list, metada_list):
                self.add_metadata_to_mol(mol, md_dict)

        with open(results_sdf_fname, "w") as f:
            mol_list = list( mol_list )
            mol_example = mol_list[0]
            missing_properties = self._autodetect_missing_scores_in_header( mol_example)
            missing_properties_names = list(map(lambda x: x[0], missing_properties))
            w = Chem.SDWriter(f)
            # create dummy mol as file header
            w.write( self._get_header_mol(missing_properties) )

            valid_props = set(list(FragalysisFormater.REQUIRED_FIELDS) + self.optional_fileds + list(missing_properties_names))


            for mol in mol_list:
                mol = load_mol_if_str(mol)
                if mol is None: continue
                mol_name = mol.GetProp("_Name")
                if not mol.HasProp("ref_pdb"):
                  assert self.ref_pdb_xchemId is not None, "Error, if ref_pdb not included in molecules as prop, self.ref_pdb_xchemId cannot be None"
                  mol.SetProp("ref_pdb", self.ref_pdb_xchemId)  # TODO: allow for more flexible future usage. see scoring_config.py

                if not mol.HasProp("original SMILES"):
                    mol.SetProp("original SMILES",  Chem.MolToSmiles(mol))

                assert self._check_molecule_metadata(mol, missing_properties_names), \
                       "Error, molecule (%s) does not contain (%s) required information (%s)." % (mol_name, str(mol.GetPropsAsDict()),
                                     str( FragalysisFormater.REQUIRED_FIELDS))

                fragments =  mol.GetProp("ref_mols")
                if self.drop_unknown_fields:
                    allProps = mol.GetPropsAsDict()
                    for prop,val in allProps.items():
                        if not prop in valid_props:
                            mol.ClearProp(prop)

                if len( fragments.strip() ) == 0:
                    journal.warning("Molecule %s has not associated fragments. Skipping..."%mol_name)
                    continue

                mol = self.trimCharFieldInMol(mol, "_Name")
                w.write(mol)
            w.close()

    def write_jsons_to_sdf(self, outname, jsons_dir, json_pattern, mols_pattern, molecules_dir=None):

        if molecules_dir is None:
            molecules_dir= jsons_dir

        mol_dict = dict( load_files_as_mols(molecules_dir, file_pattern=mols_pattern) )
        assert len( mol_dict)>0, "Error, no molecules found, please review your mols_pattern "

        def process_jsons(fname):
            with open(fname) as f:
                data = json.load(f)
            if data:
                mol_name = data["name"]
                mol = mol_dict.get(mol_name, None)
                if mol:
                    return mol, data
                return None
            else:
                return None

        per_json_result = apply_func_to_files(jsons_dir, json_pattern, process_jsons)
        assert len(per_json_result)>0, "Error, no metadata jsons found, please review your json_pattern"

        per_json_result = filter(None.__ne__, per_json_result)
        mols_to_write = (  self.add_metadata_to_mol(mol, metadata_dict) for mol, metadata_dict in per_json_result )
        self.write_molsList_to_sdf(outname, mols_to_write)
        return mols_to_write

    def write_sdfs_to_sdf(self, outname, sdfs_dir, sdfs_pattern=".sdf$"):

        sdf_names = apply_func_to_files(sdfs_dir, sdfs_pattern, lambda x: x)
        mols_to_write = []
        for fname in sdf_names:
            mols_to_write += [ mol for mol in  Chem.SDMolSupplier(fname) ]

        self.write_molsList_to_sdf(outname, mols_to_write)
        return mols_to_write

    def metadata_example(self):
        print("\n############## Metadata fields file Example ############")
        for key, val in self.metadata_header.items():
            print('"%s","%s"'%(key,val))
        print("###############################################\n")



    def filter_sdf_file(self, outname, fname_in, list_of_filters):
        '''

        :param outname:
        :param fname_in:
        :param list_of_expresions: List of tuples ("feature_name", operator_name, value)
        :return:
        '''

        assert len(list_of_filters) >0, "Error, at least one filter required"
        import operator
        ops_map = {">": operator.gt, "<": operator.lt,">=": operator.ge, "<=": operator.le}
        lambdas_list =[]
        for filt in list_of_filters:
            match_obj = re.match("(\w+)\s*(>=|<=|>|<)\s*(\d*\.?\d*)", filt)
            if match_obj:
                name, op, val = match_obj.groups()
                op = ops_map[op]
                val = float(val)
                fun = lambda mol: op(float(mol.GetProp(name)),  val)
                lambdas_list.append ( fun )
            else:
                raise  ValueError("Error, %s could not be parsed"%filt)

        do_keep = lambda mol: reduce(lambda prev, fun: prev and fun(mol), lambdas_list, True)
        mols_to_write = []

        for mol in Chem.SDMolSupplier(fname_in):
            if mol.GetProp("_Name") == FragalysisFormater.FRAGALYSIS_UPLOAD_VERSION or do_keep(mol):
                mols_to_write.append( mol)

        if  len(mols_to_write) <2:
            raise ValueError("Error, no molecules matching filters were found" )

        if outname:
            self.write_molsList_to_sdf(outname, mols_to_write)
        return mols_to_write

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser( description="Operate results to make them compatible with Fragalysis")
    subparsers = parser.add_subparsers(help='Differnt options to perform', dest='option', required=True)

    parser_merge = subparsers.add_parser('combine_sdf', help='Concatenates a set of sdf files already prepared for Fragalysis')
    parser_merge.add_argument('input_files', nargs="+" , type=str, help='sdf files already formated with Fragalysis header/labels')
    parser_merge.add_argument('-o', '--output_file', nargs=None , type=str, help='output sdf filename where the combined file will be save', required=True)

    parser_merge = subparsers.add_parser('metadata_example', help='Prints an example header')
    
    
    parser_merge = subparsers.add_parser('gather_scores', help='Creates an sdf file compatible with Fragalysis by recovring the scores contained in'
                                                               'other sdf files or .mol and .json files. It also explores subdirectories.')
    parser_merge.add_argument('data_dir', nargs=None , type=str, help='The directory where files are contained (can be located in subdirectories). '
                                                                      'It must contain either json and mol files with same name pattern or sdf files '
                                                                      'with atributes')
    parser_merge.add_argument('-o', '--output_file', nargs=None , type=str, help='output sdf filename where the combined scores will be save', required=True)
    parser_merge.add_argument('-r', '--ref_pdb_xchemId', nargs=None , type=str, help='Reference atomic model using xchem code', required=True)

    parser_merge.add_argument('-f', '--input_format',choices=['json', 'sdf'], help='The format of the input files', required=True)


    json_arg = parser_merge.add_argument('-j', '--json_pattern', nargs=None , type=str, help='regex pattern for score jsons. Only required when format is json',   required=False, default=None)
                              #,default="([\w-]+)\.scores\.json$")


    mol_arg = parser_merge.add_argument('-p', '--mols_pattern', nargs=None , type=str, help='regex pattern for .mol or .sdf files. It must contain () to match the compound id E.g. '
                                                                                  '"([\w-]+)\.minimised.mol$" ', required=False, default=None)
                              # ,default="([\w-]+)\.minimised\.mol$")


    parser_merge = subparsers.add_parser('filter', help='Filters an sdf file for Fragalysis to select a subset of molecules matching certain properties')
    parser_merge.add_argument('input_file', nargs=None , type=str, help='The sdf file to filter. It must follow Fralaysis v1.2 version')
    parser_merge.add_argument('-o', '--output_file', nargs=None , type=str, help='output sdf filename where the preserved compounds will be save', required=True)

    parser_merge.add_argument('-f', '--filters', nargs="+", type=str, help="The filers to apply in the form: 'field_name operatator' value", required=True)


    args = parser.parse_args()

    if args.option == "combine_sdf":
        FragalysisFormater.merge_sdfs(args.output_file, *args.input_files)
        '''
python -m fragmenstein.external.uploadToFragalysis.fragalysisFormater combine_sdf f1.sdf f2.sdf -o out.sdf
        '''
    elif args.option == "gather_scores":

        if args.input_format == "json":
            if args.json_pattern is None:
                raise argparse.ArgumentError(json_arg, "Error, if input_format==json, json_pattern must be provided. E.g 'scores.json$'")
            if args.mols_pattern is None:
                raise argparse.ArgumentError(mol_arg, "Error, if input_format==json, mols_pattern must be provided. E.g '(\w+).mol'")

            FragalysisFormater(args.ref_pdb_xchemId).write_jsons_to_sdf(outname=args.output_file,
                                                                        jsons_dir=args.data_dir,
                                                                        json_pattern=args.json_pattern,
                                                                        mols_pattern=args.mols_pattern)
            '''
python -m fragmenstein.external.uploadToFragalysis.fragalysisFormater gather_scores output -f json -r x0020 -o out.sdf -j '([\w-]+)\.scores\.json$' -p '([\w-]+)\.minimised\.mol$'
            '''

        elif args.input_format == "sdf":
            if args.json_pattern is not None:
                raise argparse.ArgumentError(json_arg, "Error, if input_format==sdf, json_pattern must not be provided")
            if args.mols_pattern == None:
                args.mols_pattern="(\w+)\.sdf$"

            FragalysisFormater(args.ref_pdb_xchemId).write_sdfs_to_sdf(outname=args.output_file,
                                                                       sdfs_dir=args.data_dir,
                                                                       sdfs_pattern=args.mols_pattern)
            '''
python -m fragmenstein.external.uploadToFragalysis.fragalysisFormater gather_scores sdf_trials/ -f sdf -r x0020 -o out.sdf -p '.*\.sdf$'
            '''

    elif args.option == "metadata_example":
        FragalysisFormater(None).metadata_example()
        '''
python -m fragmenstein.external.uploadToFragalysis.fragalysisFormater metadata_example
        '''

    elif args.option == "filter":
        FragalysisFormater(None).filter_sdf_file( outname=args.output_file, fname_in= args.input_file, list_of_filters=args.filters)
        '''
python -m fragmenstein.external.uploadToFragalysis.fragalysisFormater filter results_prueba.sdf -o results_prueba_filtered.sdf -f "xcos_score>0.1"

        '''
    else:
        raise ValueError("Option not valid: %s"%(args.option))
