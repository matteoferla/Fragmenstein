import os

from rdkit import Chem

import numpy as np
from Bio.PDB.PDBParser import PDBParser

class PdbDistanceManager():

    def __init__(self, pdb_fname, limit_to_resname=None):
        self.pdb_fname = pdb_fname
        self.struct = PDBParser(QUIET=True).get_structure(os.path.basename(pdb_fname), pdb_fname)
        self.pdb_residues = list(self.struct.get_residues())
        if limit_to_resname:
            self.pdb_residues = list( filter(lambda r: r.resname == limit_to_resname, self.pdb_residues))

        def getCoord(res):
            try:
                return res["CA"].coord
            except KeyError:
                try:
                    return res["CB"].coord
                except KeyError:
                    return list(res.get_atoms())[0].coord

        self.pdb_coords = np.array( [ getCoord(r) for r in self.pdb_residues] )

    def find_closest_residue(self, mol):
        center = np.mean (mol.GetConformer().GetPositions(), axis=0)
        dist2 = np.sum( (self.pdb_coords - center)**2, axis=1)
        idx = np.argmin(dist2)
        residue = self.pdb_residues[idx]
        __, __, chainId, resId = residue.get_full_id()
        resId = list(resId)
        resId[1] = str( resId[1])
        resId= "".join([ elem.strip() for elem in resId])
        return chainId, resId, residue.resname



def test():
    data_dir =  os.path.abspath( os.path.join(__file__, "../../mpro/data") )
    pdb_fname = os.path.join(data_dir, "template.pdb")
    mol =  os.path.join(data_dir, "hit_mols", "Mpro-x0678.mol")
    mol = Chem.MolFromMolFile(mol)
    found = PdbDistanceManager(pdb_fname, limit_to_resname="CYS").find_closest_residue(mol)
    print( found )

# test()