import json
import pickle
from typing import Union, OrderedDict

import logging
from rdkit import Chem
from rdkit.Chem import PropertyMol

from fragmenstein.scoring.scorer_labels import checkIfNameIsScore
from fragmenstein.utils.io_utils import load_mol

journal = logging.getLogger("Compound")

class Compound(Chem.Mol):

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

    @classmethod
    def DeleteSubstructs(cls, compound, *args, **kwargs ):
        new_mol = Chem.DeleteSubstructs( Chem.Mol(compound), *args, **kwargs)
        return cls.copyMetadata( new_mol, compound)

    @classmethod
    def GetMolFrags(cls, compound, *args, **kwargs ):
        new_mols = Chem.GetMolFrags( Chem.Mol(compound), asMols=True, *args, **kwargs)
        return [ cls.copyMetadata( new_mol, compound)  for new_mol in new_mols]

    @classmethod
    def SanitizeMol(cls, compound, *args, **kwargs ):
        new_mol = Chem.Mol(compound.ToBinary())
        Chem.SanitizeMol( new_mol, *args, **kwargs)
        return  cls.copyMetadata( new_mol, compound)

    @classmethod
    def copyMetadata(cls, mol, compound):
        comp = Compound(mol)
        comp.parents = compound.parents
        comp.molId = compound.molId
        comp.covalent_info = compound.covalent_info
        comp.scores_dict = compound.scores_dict
        comp.ref_pdb = compound.ref_pdb
        comp.ref_molIds = compound.ref_molIds
        comp.metadata = compound.metadata
        comp.unminimized_mol_pdbblock = compound.unminimized_mol_pdbblock
        return comp

    def __init__(self, mol: Union[Chem.Mol, Chem.PropertyMol.PropertyMol, "Compound"], molId=None, parents=None,
                 ref_pdb=None, *args, **kwargs):


        super().__init__(mol, *args, **kwargs)

        covalent_info = None
        metadata = None
        if isinstance(mol, type(self)):
            parents = mol.parents
            molId = molId if molId else mol.molId
            metadata = mol.metadata

        self.covalent_info = covalent_info  #E.g. {'covalent_resi':'145A', 'covalent_resn':'CYS'}
        self.props_as_dict = {}
        if parents:
            self.parents = parents
        else:
            self._parents = []
        if molId:
            self.molId = molId
        else:
            self._molId = None

        if ref_pdb:
            self.ref_pdb = ref_pdb
        else:
            self._ref_pdb = None

        self._ref_molsIds = None
        self._metadata = None

        self._scores_dict = {}

        self.disregarded_frags= None

        if metadata:
            self.metadata = metadata

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
            try:
                self._ref_pdb = self.GetProp("ref_pdb")
            except KeyError:
                return None
        return self._ref_pdb

    @ref_pdb.setter
    def ref_pdb(self, new_pdb):
        if new_pdb is None:
            journal.warning("ref_pdb is None when setting  for %s"%self)
            return
        self.SetProp("ref_pdb", new_pdb)
        self._ref_pdb = new_pdb


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
        if frag_ids is None:
            journal.warning("Warning. frag_ids is None when setting  for %s"%self)
            return
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

    @property
    def unminimized_mol_pdbblock(self):
        try:
            return self.GetProp("unminimized_mol_pdbblock")
        except KeyError:
            return None

    @unminimized_mol_pdbblock.setter
    def unminimized_mol_pdbblock(self, pdbblock):
        self.SetProp("unminimized_mol_pdbblock", pdbblock )


    def _set_typedScore(self, key, val ):
        if isinstance(val, int ):
            self.SetIntProp(key, val)
        elif isinstance(val, float):
            self.SetDoubleProp(key, val)
        else:
            self.SetProp(key, val)

    def add_scores(self, new_scores_dict):
        if new_scores_dict is None:
            journal.warning("Warning. new_scores_dict is None when setting  for %s"%self)
            return
        current_scores = self.scores_dict.copy()
        for key, val in new_scores_dict.items():
            if key in self.scores_dict:
                raise ValueError("Error, key (%s) already included  in scores"%(key))
            else:
                if checkIfNameIsScore(key):
                    self._set_typedScore(key, val)
        self.scores_dict = current_scores
        # self.SetProp("scores_dict", json.dumps(self.scores_dict ) )


    def getFragIds(self):
        return self.ref_molIds

    def __str__(self):
        return "<"+str(type(self).__name__)+": "+self.molId+">"


    __getstate_manages_dict__ = True

    def __getstate__(self):
        self.props_as_dict = self.GetPropsAsDict()
        self.mol_as_bin = self.ToBinary()

        return self.__dict__

    def __setstate__(self, stateD):
        # print( stateD )
        type(self)( stateD["mol_as_bin"] )

        self.__dict__ = stateD

        self.molId = self.molId #Required to activate SetProp("_Name", name)

        for key, val in self.props_as_dict.items():
            self._set_typedScore(key, val)

    ### This also works
    # def __getstate__(self):
    #     pDict = {}
    #     for pn in self.GetPropNames(includePrivate=True):
    #         pDict[pn] = self.GetProp(pn)
    #     return {'pkl': self.ToBinary(), 'propD': pDict, 'parents': pickle.dumps(self.parents), "molId":self.molId}
    #
    # def __setstate__(self, stateD):
    #     Chem.Mol.__init__(self, stateD['pkl'])
    #     for prop, val in stateD['propD'].items():
    #         self._set_typedScore(prop, val)
    #     self.parents = pickle.loads(stateD["parents"])
    #     self.molId = stateD["molId"]

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
    m1.molId ="kk"

    m1.parents = [m2, m3]
    print(m1.parents)

    m_str = pickle._dumps( m1)
    m1 = pickle.loads( m_str )
    print( m1.GetProp("_Name") )
    print( type(m1))
    print( m1.molId )
    print("parents after pickle", m1.parents)
    print( m1.primitiveId )

def test_pickleScores():
    m1,m2,m3,m4 = test_getMols()
    m1.add_scores({"s1_score":1})
    print( m1.scores_dict  )
    m_str = pickle.dumps( m1)
    m1 = pickle.loads( m_str )
    print( m1.scores_dict )

if __name__=="__main__":
    # test_primitiveId()
    test_pickleProps()
    # test_pickleScores()

    '''

python -m fragmenstein.protocols.dataModel.compound

    '''
