import os
import re
import zipfile

import numpy as np
from itertools import cycle

from plip.visualization.pymol import PyMOLVisualizer
from plip.structure.preparation import PDBComplex
from plip.basic.remote import VisualizerData

from rdkit import Chem
from rdkit.Chem import PandasTools

from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.protocols.steps.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.protocols.xchem_info import Xchem_info
from fragmenstein.scoring.scorer_labels import checkIfNameIsScore, removeScoreTag, SCORE_NAME_TEMPLATE
from fragmenstein.utils.io_utils import apply_func_to_files


class InteractiveInterface():

    MOL_FIELD_ID = "mol"


    EXIT_OPTION = 0
    NEXT_COMPOUND_OPTION = 1
    PREV_COMPOUND_OPTION = 2
    GO_TO_MOL_OPTION = 3
    SORT_BY_OPTION = 4
    COMBINE_SCORES = 5

    def __init__(self, sdf_fname, fragments_dir, results_dir, fragment_id_pattern=None, unboundPdb_id_pattern=None,
                 predicted_boundPdb_id_pattern=None, filter_less_than=None, filter_greater_than= None):



        with open(sdf_fname, "rb") as f:

            first_line = f.readline().decode("UTF-8")
            if first_line.startswith(FragalysisFormater.FRAGALYSIS_UPLOAD_VERSION):
                for line in f:
                    if line.decode("UTF-8").startswith("$$$$"):
                        break
            else:
                f.seek(0)

            self.df = PandasTools.LoadSDF(f, idName="molId", smilesName='smi', molColName=self.MOL_FIELD_ID, embedProps=True,
                                     includeFingerprints=False)

            self.scores_names = sorted([x for x in self.df.columns if checkIfNameIsScore(x)])
            for score_name in self.scores_names:
                self.df[score_name] = self.df[score_name].astype(np.float32)

        if filter_less_than:
            for key_val in filter_less_than:
                key, val = key_val.split(":")
                key = key+"_score"
                self.df = self.df[ self.df[key] < float(val)]

        if filter_greater_than:
            for key_val in filter_greater_than:
                key, val = key_val.split(":")
                key = key+"_score"
                self.df = self.df[ self.df[key] > float(val)]
        self.df.reset_index(inplace=True)

        # self.df = self.df.sort_values(by=self.scores_names[0], ascending=True)

        if not fragment_id_pattern:
            fragment_id_pattern = LoadInput_XchemDefault.fragment_id_pattern
        if not unboundPdb_id_pattern:
            unboundPdb_id_pattern = LoadInput_XchemDefault.unboundPdb_id_pattern
        if not predicted_boundPdb_id_pattern:
            predicted_boundPdb_id_pattern = LoadInput_XchemDefault.predicted_boundPdb_id_pattern

        self.data_loader = LoadInput_XchemDefault( fragments_dir, fragment_id_pattern=fragment_id_pattern,
                                                  unboundPdb_id_pattern=unboundPdb_id_pattern)

        self.minmized_templates ={}
        with zipfile.ZipFile(os.path.join(results_dir, "bound_minimized_pdbs.zip")) as zip_f:
            # with zip_f.open('bound_pdbs/compoundName_to_pdbName.tab') as text_f:
            with zip_f.open('info/compoundName_to_pdbName.tab') as text_f:
                for line in text_f:
                    line = line.decode("utf-8").split()
                    if len(line)<2:
                        continue
                    self.minmized_templates[line[0]] = os.path.join(results_dir, "merges", line[-1],
                                                                    Xchem_info.predicted_boundPdb_template%line[-1])

        self.molecule_idx=0
        self.n_mols = len(self.df)

        cls = type(self)
        self.options_dict = {
            r"exit": (cls.EXIT_OPTION, lambda groups: True), # sys.exit(0)),
            r"(next|\n)": (cls.NEXT_COMPOUND_OPTION, lambda groups: self.next_mol() ),
            r"prev": (cls.NEXT_COMPOUND_OPTION, lambda groups: self.prev_mol()),
            r"goto (-?\d+)$": (cls.GO_TO_MOL_OPTION, lambda groups: self.goto(groups)),
            r"jumpto ([-_\w]+)$": (cls.GO_TO_MOL_OPTION, lambda groups: self.jump(groups)),
            r"sortby (\w+) (asc|desc)": (cls.SORT_BY_OPTION, lambda groups: self.sortby(groups)),
            r"combinescores (\w+) (.+)": (cls.COMBINE_SCORES, lambda groups: self.combineScores(groups)),

        }

    def next_mol(self, *args):

        if self.molecule_idx >= (self.n_mols-1):
            print("No more molecules to display! Current: %d Total: %d"%(self.molecule_idx+1, self.n_mols))
        else:
            self.molecule_idx+= 1

    def prev_mol(self, *args):

        if self.molecule_idx <=0:
            print("No more molecules to display! Current: %d Total: %d"%(self.molecule_idx+1, self.n_mols))
        else:
            self.molecule_idx-= 1

    def goto(self, num_match):

        num = int( num_match[0])
        if num<0:
            new_idx = self.n_mols + num
        else:
            new_idx = num -1

        if new_idx <0 or new_idx > (self.n_mols-1):
            print("No more molecules to display! Current: %d Total: %d"%(self.molecule_idx+1, self.n_mols))
        else:
            self.molecule_idx = new_idx

    def sortby(self, feature_order):

        feature, order = feature_order
        if feature not in self.scores_names:
            feature = SCORE_NAME_TEMPLATE%feature
            if feature not in self.scores_names:
                print("Error, selected feature does not match known features: %s" % str(self.scores_names))
                return

        self.df = self.df.sort_values(by=feature, ascending= True if order.startswith("asc") else False)
        self.molecule_idx = 0

    def combineScores(self, name_expresion):
        '''

        E.g.         'combinescores s1 plipGlobalPreser_score/rotableBonds_score'

        :param name_expresion:
        :return:
        '''
        name, expression = name_expresion
        if not checkIfNameIsScore(name):
            name = SCORE_NAME_TEMPLATE%name
        if name in self.df:
            print("Error, feature name already exists: %s" % str(name))
            return
        self.df[name] = self.df.eval(expression)
        def updateMol(row):
            newFeat = row[name]
            mol = row[self.MOL_FIELD_ID]
            mol.SetProp(name, str(newFeat))
            return mol

        self.df[ self.MOL_FIELD_ID] = self.df.apply( updateMol, axis=1)
        self.scores_names.append( name)
        pass

    def get_current_mol(self):
        name, mol = self.df.iloc[self.molecule_idx, :][["molId", "mol"]]
        return name, mol

    def get_fragments_current_mol(self):
        frag_ids = self.get_current_mol()[-1].GetProp("ref_mols").split(",")
        return [(fragId, self.data_loader.fragments_dict[fragId]) for fragId in frag_ids]


    def get_mol_repr(self, mol):

        # print( mol.ToBinary())

        sorted_tuples = sorted([(removeScoreTag(key), round(val, 6)) for key, val in mol.GetPropsAsDict().items() if
                                checkIfNameIsScore(key)], key=lambda x: len(x[0]))
        score_str = mol.GetProp("_Name")+" -> ["+mol.GetProp("ref_mols")+"]\n"
        score_str += Chem.MolToSmiles(mol)+"\n"
        template = "%10s: %9s\t"
        for name, val in sorted_tuples:
            score_str += template % (name, val)
        return score_str

    def ask_option(self):
        print("#######################################")
        print("# Current: %d Total: %d"%(self.molecule_idx+1, self.n_mols))
        print("# Options are:")
        for option in self.options_dict:
            print("#   - "+repr(option))
        print("#######################################")
        try:
            users_in = input("Type the selected option: ",)
            print(".......................................")
            if users_in=="":
                users_in = "\n"
            else:
                users_in = users_in.strip()

            for option in self.options_dict:
                match_obj = re.match(option, users_in)
                if match_obj:
                    return self.options_dict[option][-1](match_obj.groups())
        except UnicodeDecodeError:
            pass
        print("***************************************")
        print("Incorrect option!")
        return self.ask_option()

    def jump(self, groups):
        molId = groups[0]
        indices = self.df.index[self.df[ "molId" ] == molId].tolist()
        if not indices:
            print("Error id (%s) not found! Current: %d Total: %d" % (molId, self.molecule_idx + 1, self.n_mols))
        else:
            self.molecule_idx = indices[0]

    def display_interactions(self, pdbfname):

        example = PDBComplex()
        example.load_pdb(pdbfname)
        example.analyze()
        plcomplex = None
        for site in sorted(example.interaction_sets):
            if not len(example.interaction_sets[site].interacting_res) == 0 and site.startswith("LIG"):
                plcomplex = VisualizerData(example, site)
                break
        if not plcomplex:
            return

        vis = PyMOLVisualizer(plcomplex)

        #####################
        # Set everything up #
        #####################

        pdbid = plcomplex.pdbid
        lig_members = plcomplex.lig_members
        chain = plcomplex.chain
        from plip.basic import config
        if config.PEPTIDES:
            vis.ligname = 'PeptideChain%s' % plcomplex.chain
        if config.INTRA is not None:
            vis.ligname = 'Intra%s' % plcomplex.chain

        ligname = vis.ligname
        hetid = plcomplex.hetid

        metal_ids = plcomplex.metal_ids
        metal_ids_str = '+'.join([str(i) for i in metal_ids])

        ########################
        # Basic visualizations #
        ########################
        from pymol import cmd

        vis.set_initial_representations()

        cmd.load(plcomplex.sourcefile)
        cmd.frame(config.MODEL)
        current_name = cmd.get_object_list(selection='(all)')[0]

        cmd.set_name(current_name, pdbid)
        cmd.hide('everything', 'all')
        if config.PEPTIDES:
            cmd.select(ligname, 'chain %s and not resn HOH' % plcomplex.chain)
        else:
            cmd.select(ligname, 'resn %s and chain %s and resi %s*' % (hetid, chain, plcomplex.position))

        # Visualize and color metal ions if there are any
        if not len(metal_ids) == 0:
            vis.select_by_ids(ligname, metal_ids, selection_exists=True)
            cmd.show('spheres', 'id %s and %s' % (metal_ids_str, pdbid))

        # Additionally, select all members of composite ligands
        if len(lig_members) > 1:
            for member in lig_members:
                resid, chain, resnr = member[0], member[1], str(member[2])
                cmd.select(ligname, '%s or (resn %s and chain %s and resi %s)' % (ligname, resid, chain, resnr))

        cmd.show('sticks', ligname)
        cmd.color('myblue')
        cmd.color('myorange', ligname)
        cmd.util.cnc('all')
        if not len(metal_ids) == 0:
            cmd.color('hotpink', 'id %s' % metal_ids_str)
            cmd.hide('sticks', 'id %s' % metal_ids_str)
            cmd.set('sphere_scale', 0.3, ligname)
        cmd.deselect()

        vis.make_initial_selections()

        vis.show_hydrophobic()  # Hydrophobic Contacts
        vis.show_hbonds()  # Hydrogen Bonds
        vis.show_halogen()  # Halogen Bonds
        vis.show_stacking()  # pi-Stacking Interactions
        vis.show_cationpi()  # pi-Cation Interactions
        vis.show_sbridges()  # Salt Bridges
        vis.show_wbridges()  # Water Bridges
        vis.show_metal()  # Metal Coordination

        vis.refinements()

        vis.zoom_to_ligand()

        vis.selections_cleanup()

        vis.selections_group()
        vis.additional_cleanup()
        if config.DNARECEPTOR:
            # Rename Cartoon selection to Line selection and change repr.
            cmd.set_name('%sCartoon' % plcomplex.pdbid, '%sLines' % plcomplex.pdbid)
            cmd.hide('cartoon', '%sLines' % plcomplex.pdbid)
            cmd.show('lines', '%sLines' % plcomplex.pdbid)


class PymolManager():

    templatePymolId = "template"
    compound_color = 'wheat'
    fragment_colours = cycle(
        ['palegreen', 'lightblue', 'lightpink', 'palecyan', 'bluewhite'])

    def __init__(self):

        import pymol
        self.pymol = pymol
        self.pymol.finish_launching(['pymol', '-q'])
        self.active_molIds = []

    def __del__(self): #TODO: change it using atexit

        self.pymol.cmd.quit()
        del self.pymol

    def load_atomic_model(self, fname, pymolTemplateId=None):
        if pymolTemplateId is None:
            pymolTemplateId = PymolManager.templatePymolId

        self.pymol.cmd.delete( pymolTemplateId)

        self.pymol.cmd.load(fname, pymolTemplateId)
        self.pymol.cmd.select("LIG", "resname LIG")
        self.pymol.cmd.delete("LIG")
        self.pymol.cmd.show_as("cartoon", pymolTemplateId)
        self.pymol.cmd.show("lines", pymolTemplateId)


    def draw_mol(self, name, mol, is_fragment=False):

        self.pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), name)
        if is_fragment:
            color = next(PymolManager.fragment_colours)
        else:
            color = PymolManager.compound_color

        self.pymol.cmd.color(color, '{name} and name C*'.format(name=name))

        if is_fragment:
            self.pymol.cmd.set("stick_transparency", "0.5", name)
            self.pymol.cmd.set_bond("stick_radius", "0.2", name)
            self.pymol.cmd.zoom(name, 5)
        else:
            self.pymol.cmd.select(name)
            # self.pymol.cmd.do(" sele template") # contacts template, sele")

        self.active_molIds.append( name )

    def clean_mols(self, all=True):
        if all:
            self.pymol.cmd.delete("all")
        else:
            for molId in self.active_molIds:
                self.pymol.cmd.delete(molId)


# def main(sdf_fname, template, fragments_dir, results_dir, *args, **kwargs):
#
#     sdf_fname = os.path.expanduser(sdf_fname)
#     fragments_dir = os.path.expanduser(fragments_dir)
#
#     interactiv_sess = InteractiveInterface(sdf_fname, fragments_dir, results_dir)
#     pm = PymolManager()
#
#     if template:
#         template = os.path.expanduser(template)
#         pm.load_atomic_model(template)
#
#     finished =False
#     while not finished:
#         name, mol = interactiv_sess.get_current_mol()
#         print(interactiv_sess.get_mol_repr(mol))
#         pm.draw_mol(name, mol)
#         minimizedPdbName = interactiv_sess.minmized_templates[name]
#         pm.load_atomic_model(minimizedPdbName, "minimized_model")
#         names_mols_list = interactiv_sess.get_fragments_current_mol()
#         for name, mol in names_mols_list:
#             pm.draw_mol( name, mol, is_fragment=True)
#
#         interactiv_sess.display_interactions(minimizedPdbName)
#
#         finished = interactiv_sess.ask_option()
#         # pm.clean_mols()
#         pm.pymol.cmd.delete("all")
#
#     print("Closing pymol")
#     del pm

def main(sdf_fname, template, fragments_dir, results_dir, filter_less_than=None, filter_greater_than= None, *args, **kwargs):

    sdf_fname = os.path.expanduser(sdf_fname)
    fragments_dir = os.path.expanduser(fragments_dir)

    interactiv_sess = InteractiveInterface(sdf_fname, fragments_dir, results_dir, filter_less_than=filter_less_than,
                                           filter_greater_than= filter_greater_than)
    pm = PymolManager()


    finished =False
    while not finished:
        name, mol = interactiv_sess.get_current_mol()
        print(interactiv_sess.get_mol_repr(mol))

        minimizedPdbName = interactiv_sess.minmized_templates[name]
        interactiv_sess.display_interactions(minimizedPdbName)

        # pm.draw_mol(name, mol)
        # pm.load_atomic_model(minimizedPdbName, "minimized_model")
        names_mols_list = interactiv_sess.get_fragments_current_mol()
        for name, mol in names_mols_list:
            pm.draw_mol( name, mol, is_fragment=True)

        if template:
            template = os.path.expanduser(template)
            pm.load_atomic_model(template)

        finished = interactiv_sess.ask_option()
        pm.clean_mols(all=True)

    print("Closing pymol")
    del pm


if __name__ == "__main__":
    from fragmenstein.utils.cmd_parser import ArgumentParser
    parser = ArgumentParser(prog="view_results_pymol", description="visualize molecules and scores using pymol")
    parser.add_argument("-i", "--sdf_fname", type=str, help="The file with annotated molecules ", required=True)
    parser.add_argument("-d", "--fragments_dir", type=str, help="The root dir where fragments are, typically target_name/aligned/ ", required=True)
    parser.add_argument("-r", "--results_dir", type=str, help="The path to computed results", required=True)
    parser.add_argument("-t", "--template", type=str, help="The path to a template pdb to display (in addition to minimized pdb", required=False)
    parser.add_argument("-l", "--filter_less_than", type=str, nargs='+', help="A list of key:value for variable constrains"
                                  " to be satisfied so that variables[key]<value. E.g. energyPoseVsEnsemble:100", required=False)
    parser.add_argument("-g", "--filter_greater_than", type=str, nargs='+', help="A list of key:value for variable constrains"
                                  " to be satisfied so that variables[key]>value. E.g. sumSuCosW:0.5", required=False)

    args =vars( parser.parse_args())
    print(args)
    main( ** args)
    print("\nmain DONE!\n")

'''

python -m fragmenstein.external.view_results_pymol.view_results_pymol -i ~/oxford/myProjects/diamondCovid/data/nsp13/enumeration/Site_1/site_1_perm.sdf -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -t ~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0020_0B/nsp13-x0020_0B_apo-desolv.pdb

'''


'''
contacts template, sele

'''