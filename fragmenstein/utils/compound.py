from typing import Union

from rdkit import Chem
from rdkit.Chem import PropertyMol

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

    def __init__(self, mol: Union[Chem.Mol, Chem.PropertyMol.PropertyMol], molId=None, parents=None, *args, **kwargs):

        super().__init__(mol, *args, **kwargs)

        self.covalent_info = None  #E.g. {'covalent_resi':'145A', 'covalent_resn':'CYS'}

        if parents:
            self._parents = parents
        else:
            self._parents = []
        if molId:
            self._molId = molId
        else:
            self._molId = None

    @property
    def parents(self):
        return self._parents

    @parents.setter
    def parents(self, parents):
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

def test_pickleProps():
    import  pickle
    m1,m2,m3,m4 = test_getMols()
    m1.SetProp("_Name", "kk")
    m_str = pickle._dumps( m1)
    m1 = pickle.loads( m_str )
    print( m1.GetProp("_Name") )
    print( type(m1))
    print( m1.molId )


if __name__=="__main__":
    test_primitiveId()
    test_pickleProps()