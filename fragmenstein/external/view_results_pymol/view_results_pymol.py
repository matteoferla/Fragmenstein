import os
import re

from rdkit.Chem import PandasTools
from itertools import cycle

from rdkit import Chem


from fragmenstein.external.uploadToFragalysis.fragalysisFormater import FragalysisFormater
from fragmenstein.protocols.loadInput_XchemDefault import LoadInput_XchemDefault
from fragmenstein.scoring.scorer_labels import checkIfNameIsScore, removeScoreTag, SCORE_NAME_TEMPLATE


class InteractiveInterface():

    EXIT_OPTION = 0
    NEXT_COMPOUND_OPTION = 1
    PREV_COMPOUND_OPTION = 2
    GO_TO_MOL_OPTION = 3
    SORT_BY_OPTION = 4

    def __init__(self, sdf_fname, fragments_dir, fragment_id_pattern=None, unboundPdb_id_pattern=None):



        with open(sdf_fname, "rb") as f:

            first_line = f.readline().decode("UTF-8")
            if first_line.startswith(FragalysisFormater.FRAGALYSIS_UPLOAD_VERSION):
                for line in f:
                    if line.decode("UTF-8").startswith("$$$$"):
                        break
            else:
                f.seek(0)

            self.df = PandasTools.LoadSDF(f, idName="molId", smilesName='smi', molColName='mol', embedProps=True,
                                     includeFingerprints=False)

            self.scores_names = sorted([x for x in self.df.columns if checkIfNameIsScore(x)])


        self.df = self.df.sort_values(by=self.scores_names[0], ascending=True)

        if not fragment_id_pattern:
            fragment_id_pattern = LoadInput_XchemDefault.fragment_id_pattern
        if not unboundPdb_id_pattern:
            unboundPdb_id_pattern = LoadInput_XchemDefault.unboundPdb_id_pattern

        self.data_loader = LoadInput_XchemDefault( fragments_dir, fragment_id_pattern=fragment_id_pattern,
                                                  unboundPdb_id_pattern=unboundPdb_id_pattern)

        self.molecule_idx=0
        self.n_mols = len(self.df)

        cls = type(self)
        self.options_dict = {
            r"exit": (cls.EXIT_OPTION, lambda groups: True), # sys.exit(0)),
            r"(next|\n)": (cls.NEXT_COMPOUND_OPTION, lambda groups: self.next_mol() ),
            r"prev": (cls.NEXT_COMPOUND_OPTION, lambda groups: self.prev_mol()),
            r"goto (-?\d+)$": (cls.GO_TO_MOL_OPTION, lambda groups: self.goto(groups)),
            r"sortby (\w+) (asc|desc)": (cls.SORT_BY_OPTION, lambda groups: self.sortby(groups)),

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

    def get_current_mol(self):
        name, mol = self.df.iloc[self.molecule_idx, :][["molId", "mol"]]
        return name, mol

    def get_fragments_current_mol(self):
        frag_ids = self.get_current_mol()[-1].GetProp("ref_mols").split(",")
        return [(fragId, self.data_loader.fragments_dict[fragId]) for fragId in frag_ids]


    def get_mol_repr(self, mol):
        sorted_tuples = sorted([(removeScoreTag(key), round(val, 6)) for key, val in mol.GetPropsAsDict().items() if
                                checkIfNameIsScore(key)], key=lambda x: len(x[0]))
        score_str = mol.GetProp("_Name")+" -> ["+mol.GetProp("ref_mols")+"]\n"
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
        print("***************************************")
        print("Incorrect option!")
        return self.ask_option()



class PymolManager():

    templatePymolId = "template"
    compound_color = 'wheat'
    fragment_colours = cycle(
        ['palegreen', 'lightblue', 'paleyellow', 'lightpink', 'palecyan', 'lightorange', 'bluewhite'])

    def __init__(self):

        import pymol
        self.pymol = pymol
        self.pymol.finish_launching(['pymol', '-q'])
        self.active_molIds = []

    def load_atomic_model(self, fname):
        self.pymol.cmd.load(fname, PymolManager.templatePymolId)
        self.pymol.cmd.show_as("cartoon", PymolManager.templatePymolId)
        self.pymol.cmd.show("sticks", PymolManager.templatePymolId)

    def __del__(self):

        self.pymol.cmd.quit()
        del self.pymol


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

        self.active_molIds.append( name )

    def clean_mols(self):
        for molId in self.active_molIds:
            self.pymol.cmd.delete(molId)


def main(sdf_fname, template, fragments_dir, *args, **kwargs):

    sdf_fname = os.path.expanduser(sdf_fname)
    fragments_dir = os.path.expanduser(fragments_dir)
    template = os.path.expanduser(template)

    interactiv_sess = InteractiveInterface(sdf_fname, fragments_dir)
    pm = PymolManager()

    pm.load_atomic_model(template)

    finished =False
    while not finished:
        name, mol = interactiv_sess.get_current_mol()
        print(interactiv_sess.get_mol_repr(mol))
        pm.draw_mol(name, mol)
        names_mols_list = interactiv_sess.get_fragments_current_mol()
        for name, mol in names_mols_list:
            pm.draw_mol( name, mol, is_fragment=True)
        finished = interactiv_sess.ask_option()
        pm.clean_mols()

    print("Closing pymol")
    del pm


if __name__ == "__main__":
    from fragmenstein.utils.cmd_parser import ArgumentParser
    parser = ArgumentParser(prog="view_results_pymol", description="visualize molecules and scores using pymol")
    parser.add_argument("-i", "--sdf_fname", type=str, help="The file with annotated molecules ", required=True)
    parser.add_argument("-d", "--fragments_dir", type=str, help="The root dir where fragments are, typically target_name/aligned/ ", required=True)
    parser.add_argument("-t", "--template", type=str, help="The path to a template pdb", required=False) #TODO: find template within minizmied.pdb

    args =vars( parser.parse_args())
    print(args)
    main( ** args)
    print("\nmain DONE!\n")

'''

python -m fragmenstein.external.view_results_pymol.view_results_pymol -i ~/oxford/myProjects/diamondCovid/data/nsp13/enumeration/Site_1/site_1_perm.sdf -d ~/oxford/myProjects/diamondCovid/data/nsp13/aligned -t ~/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0020_0B/nsp13-x0020_0B_apo-desolv.pdb

'''