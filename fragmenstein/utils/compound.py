import json
import pickle
from typing import Union, OrderedDict

from rdkit import Chem
from rdkit.Chem import PropertyMol

from fragmenstein.scoring.scorer_labels import checkIfNameIsScore
from fragmenstein.utils.io_utils import load_mol

#TODO; create module dataModel

class Compound(Chem.PropertyMol.PropertyMol):

    ID_SEP = "-"

    @classmethod
    def MolFromSmiles(cls, smi, molId=None):
        comp =  cls( Chem.MolFromSmiles( smi ))
        if molId:
            comp.molId = molId
        return comp

    @classmethod
    def MolFromFile(cls, fname, molId=None):
        comp = cls( load_mol( fname) )
        if molId:
            comp.molId = molId
        return comp

    # @classmethod
    # def DeleteSubstructs(cls, compound, *args, **kwargs ):
    #     new_mol = Chem.DeleteSubstructs(compound, *args, **kwargs)
    #     return cls.copyMetadata( new_mol, compound)

    @classmethod
    def copyMetadata(cls, mol, compound):
        comp = Compound(mol)
        comp.parents = compound.parents
        comp.molId = compound.parents
        return comp

    def __init__(self, mol: Union[Chem.Mol, Chem.PropertyMol.PropertyMol, "Compound"], molId=None, parents=None,
                 ref_pdbId=None, *args, **kwargs):


        super().__init__(mol, *args, **kwargs)

        covalent_info = None
        if isinstance(mol, type(self)):
            parents = mol.parents
            molId = molId

        self.covalent_info = covalent_info  #E.g. {'covalent_resi':'145A', 'covalent_resn':'CYS'}

        if parents:
            self.parents = parents
        else:
            self._parents = []
        if molId:
            self.molId = molId
        else:
            self._molId = None

        if ref_pdbId:
            self.ref_pdb = ref_pdbId
        else:
            self._ref_pdb = None

        self._ref_molsIds = None
        self._metadata = None

        self._scores_dict = {}

        self.disregarded_frags= None

    @property
    def parents(self):
        return self._parents

    @parents.setter
    def parents(self, parents):
        self.SetProp("parents", "-".join([elem.molId for elem in parents]))
        self._parents = parents

    @property
    def molId(self):
        if self.HasProp("_Name"):
            name = self.GetProp("_Name")
            if name != "":
                self.molId = name
        return self._molId

    @molId.setter
    def molId(self, molId):
        self.SetProp("_Name", molId)
        self._molId = molId

    @property
    def primitiveId(self):
        primary_ids = []
        if len(self.parents)>0:
            for parent in self.parents:
                primary_ids+= [ parent.primitiveId ]
            return (Compound.ID_SEP).join(sorted(set(primary_ids)))
        else:
            return self.molId

    @property
    def ref_pdb(self):
        if not self._ref_pdb:
            self._ref_pdb = self.GetProp("ref_pdb")
        return self._ref_pdb

    @ref_pdb.setter
    def ref_pdb(self, ref_pdb):
        self.SetProp("ref_pdb", ref_pdb)
        self._ref_pdb = ref_pdb


    @property
    def metadata(self):
        if not self._metadata:
            if self.HasProp("metadata"):
                self._metadata = json.loads(self.GetProp("metadata"))
        return self._metadata

    @metadata.setter
    def metadata(self, metadata_dict):
        self._metadata = metadata_dict
        self.add_scores(metadata_dict)
        self.SetProp("metadata", json.dumps( metadata_dict ) )


    @property
    def ref_molIds(self):
        if not self._ref_molsIds:
            if self.HasProp("ref_mols"):
                self._ref_molsIds = self.GetProp("ref_mols").split(",")
        return self._ref_molsIds

    @ref_molIds.setter
    def ref_molIds(self, frag_ids):
        self.SetProp("ref_mols", ",".join(frag_ids))
        self._ref_molsIds = frag_ids

    @property
    def scores_dict(self):
        if not self._scores_dict:
            if self.HasProp("scores_dict"):
                self._scores_dict = json.loads(self.GetProp("scores_dict"))
        return  self._scores_dict

    @scores_dict.setter
    def scores_dict(self, scores_dict):
        self._scores_dict = scores_dict
        self.SetProp("scores_dict", json.dumps(self._scores_dict ) )

    def add_scores(self, scores_dict):
        for key, val in scores_dict.items():
            if key in self.scores_dict:
                raise ValueError("Error, key (%s) already included  in scores"%(key))
            else:
                if checkIfNameIsScore(key):
                    self.SetProp(key, val)
                    self.scores_dict[key] = val
        self.SetProp("scores_dict", json.dumps(self.scores_dict ) )


    def getFragIds(self):
        return self.ref_molIds

    def __str__(self):
        return "<"+str(type(self).__name__)+": "+self.molId+">"

def test_getMols():
    m1 = Compound.MolFromSmiles("CC")
    m1.molId = "m1"
    m2 = Compound.MolFromSmiles("CCO")
    m2.molId = "m2"

    m3 = Compound.MolFromSmiles("CCCO")
    m3.molId = "m3"

    m4 = Compound.MolFromSmiles("CCCCO")
    m4.molId = "m4"
    return m1,m2,m3,m4

def test_primitiveId():
    m1,m2,m3,m4 = test_getMols()

    # m3.parents = [m1]
    # print( m3.parents)
    # print(m3.primitiveId)

    # m3.parents = [m1, m2]
    # print( m3.parents)
    # print(m3.primitiveId)

    # m2.parents = [m1]
    # m3.parents = [m1, m2]
    # print( m3.parents)
    # print(m3.primitiveId)

    # m4.parents = [m3, m2]
    # m3.parents = [m1 ]
    # print(m4.primitiveId)

    m4.parents = [m3]
    m3.parents = [m2]
    m2.parents = [m1]
    print(m4.primitiveId)
    Compound(m2)
def test_pickleProps():
    import  pickle
    m1,m2,m3,m4 = test_getMols()
    m1.SetProp("_Name", "kk")
    m_str = pickle._dumps( m1)
    m1 = pickle.loads( m_str )
    print( m1.GetProp("_Name") )
    print( type(m1))
    print( m1.molId )

def test_pickleScores():
    m1,m2,m3,m4 = test_getMols()
    m1.add_scores({"s1_score":1})
    print( m1.scores_dict  )
    m_str = pickle.dumps( m1)
    m1 = pickle.loads( m_str )
    print( m1.scores_dict )

if __name__=="__main__":
    # test_primitiveId()
    # test_pickleProps()
    test_pickleScores()

    '''

python -m fragmenstein.utils.compound

    '''
