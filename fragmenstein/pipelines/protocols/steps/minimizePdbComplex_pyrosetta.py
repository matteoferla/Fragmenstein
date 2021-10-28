import json
import os
import numpy as np

from tempfile import NamedTemporaryFile
from typing import List, Callable

from rdkit.Chem import rdFMCS
from rdkit_to_params import constraint, Params

from rdkit import Chem
from rdkit.Chem import AllChem

from fragmenstein import Igor

from fragmenstein.external import ExternalToolImporter
from fragmenstein.monster import GPM
from fragmenstein.utils.config_manager import ConfigManager
from fragmenstein.utils.pdb_utils import PdbDistanceManager
from fragmenstein.utils.timeout import timeout
from fragmenstein.victor import MinimalPDBParser


class MinimizePDBComplex_pyrosetta():


    LIGAND_PLACED_RESID= (1, "Z")
    LIGAND_PLACED_RESNAME = "LIG"

    def __init__(self, templateFname, atom_constrain_filter, restrain_type='HARMONIC 0 1',
        *args, **kwargs):

        super().__init__( *args, **kwargs)

        self.atom_constrain_filter = atom_constrain_filter # e.g.   lambda atom: atom.HasProp("is_original_atom")
        self.templateFname = templateFname
        self.restrain_type = restrain_type

        self._pyrosetta = None
        self.igor = None

    def _appendParametrizedMolToPdbBlock(self, params, unboundFname, ): #TODO: move it to pdb_utils and make it more rosbust

        '''

        TODO: https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/DAB4155708F8E046A9D670513900A98C23200EA2%40INHEXMB09.eu.boehringer.com/#msg36404236 shows how to create residue information

        :param mol:
        :param unboundFname:
        :return:
        '''

        #TODO: DEAL WITH CL, MG and other ions.

        mol = params.mol

        lig_redId_chainId = type(self).LIGAND_PLACED_RESID
        lig_resname = type(self).LIGAND_PLACED_RESNAME

        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            info.SetResidueNumber(lig_redId_chainId[0])
            info.SetChainId( lig_redId_chainId[1] )
            info.SetIsHeteroAtom(True)
            info.SetOccupancy(1.)
            info.SetResidueName(lig_resname)

        with open(unboundFname) as f:
            pdbdata = MinimalPDBParser(f.read())

        mol = Chem.Mol( mol.ToBinary() )
        moldata = MinimalPDBParser(Chem.MolToPDBBlock(mol))
        pdbdata.append(moldata)
        return str(pdbdata)

    def get_constrains(self, params: Params, atom_constrain_filter:Callable, reference_chainId_redId_resname=None):
        '''

        :param params:
        :param atom_constrain_filter:  e.g. lambda atom: atom.HasProp("is_original_atom")
        :param reference_chainId_redId_resname:
        :return:
        '''
        mol = params.mol
        if not reference_chainId_redId_resname:
            reference_chainId_redId_resname = PdbDistanceManager(self.templateFname).find_closest_residue(mol)
        chainId, resId, resname = reference_chainId_redId_resname
        ligand_resi = str(self.LIGAND_PLACED_RESID[0])+self.LIGAND_PLACED_RESID[1]
        covalent_resi = resId + chainId
        constrains = constraint.Constraints.mock()

        conf = mol.GetConformer()
        lines = []
        n_constrained_atoms = 0
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            if atom_constrain_filter(atom):
                if atom.GetSymbol() == '*':
                    continue
                elif atom.GetPDBResidueInfo() is None:
                    print('Atom {i} ({sym}) has no name!'.format(i=i,sym=atom.GetSymbol()))
                    continue
                n_constrained_atoms += 1
                pos = conf.GetAtomPosition(i)
                atomname = atom.GetPDBResidueInfo().GetName()
                lines.append(('CoordinateConstraint {atomname} {ligand_resi} '.format(atomname=atomname, ligand_resi=ligand_resi) +
                             'CA {covalent_resi} '.format(covalent_resi=covalent_resi) +
                             " ".join(map(str,pos)) + ' {fxn}\n'.format(fxn=self.restrain_type)))

        constrains.custom_constraint += ''.join(lines)

        n_unconstrained_atoms = mol.GetNumAtoms() - n_constrained_atoms
        return  constrains, (n_constrained_atoms, n_unconstrained_atoms)


    def minimize_recipe(self, params):
        dG_unbound = self.igor.detailed_scores(params.test(), 1)['total_score']
        return self.minimize_victor_igor(dG_unbound)


    def minimize_victor_igor(self, dG_unbound) -> float:
        """
        Calls Igor recursively until the ddG is negative or zero.
        igor.minimise does a good job. this is just to get everything as a normal molecule

        :return: ddG (kcal/mol)
        """

        ddG = 999
        self.igor.coordinate_constraint = 0.
        # self.igor.fa_intra_rep = 0.02 # 4x
        # quick unconstrained minimisation to wiggle it out of nasty local minima
        self.igor.minimise(cycles=15, default_coord_constraint=False)
        self.igor.coordinate_constraint = 2
        self.igor.minimise(cycles=5, default_coord_constraint=False)
        self.igor.coordinate_constraint = 1
        while ddG > 0:
            self.igor.minimise(default_coord_constraint=False)
            dG_bound = self.igor.ligand_score()['ligand_ref2015']['total_score']
            ddG = dG_bound - dG_unbound
            if ddG > 0:
                self.igor.coordinate_constraint /= 2
            if self.igor.coordinate_constraint == 0.:
                break
            elif self.igor.coordinate_constraint < 0.005:
                self.igor.coordinate_constraint = 0.

        return ddG

    def simple_mRSMD(self, minimized_mol, placed_mol, fragments):
        '''
        :param mol:
        :param fragments:
        :param atom_constrain_filter:
        :return:
        '''

        coors = minimized_mol.GetConformer().GetPositions()

        sdevs = []
        for frag in fragments:
            mapping = GPM.get_positional_mapping(placed_mol, frag)
            atoms_idxs_mol, atoms_idxs_frag = zip(* mapping.items() )
            dev =  np.sum( (coors[atoms_idxs_mol, :] - frag.GetConformer().GetPositions()[atoms_idxs_frag, :])**2, axis=-1)
            dev = np.mean( dev)
            sdevs.append( dev )

        return np.sqrt( np.mean( sdevs ) )

    def minimize(self, mol, molId=None, outdir=None, reference_fragments=None): #TODO: clean all this thing


        if not self._pyrosetta:
            self._pyrosetta = ExternalToolImporter.import_tool("pyrosetta", ["pyrosetta"])

        Params.log.setLevel("ERROR")
        params = Params.from_mol( mol, name= type(self).LIGAND_PLACED_RESNAME)
        constrains, (n_constrained_atoms, n_unconstrained_atoms) = self.get_constrains(params, atom_constrain_filter=self.atom_constrain_filter)

        boundPdbStr = self._appendParametrizedMolToPdbBlock(params, self.templateFname )

        with NamedTemporaryFile() as tmp1, NamedTemporaryFile() as tmp2:

            params_filename = tmp1.name
            constrains_filename = tmp2.name

            params.dump(params_filename)
            constrains.dump(constrains_filename)
            tmp1.seek(0)
            tmp2.seek(0)

            self.igor = Igor.from_pdbblock(
                                pdbblock=boundPdbStr, params_file= params_filename ,
                                constraint_file= constrains_filename, ligand_residue= type(self).LIGAND_PLACED_RESNAME
            )

            unminimizedPdbBlock = self.igor.pose2str()

            ddG = self.minimize_recipe( params)


            ligand = self.igor.mol_from_pose()
            template = AllChem.DeleteSubstructs(mol, Chem.MolFromSmiles('*'))
            AllChem.AssignBondOrdersFromTemplate(template, ligand)



            if reference_fragments:
                rmsd = self.simple_mRSMD(ligand, mol, reference_fragments)

            metadata = { "∆∆G": ddG, "N_constrained_atoms": n_constrained_atoms, "N_unconstrained_atoms": n_unconstrained_atoms,
                         "comRMSD": rmsd}
            # input( metadata )
            # print(self.igor.ligand_score())
            # print("minimization done!!!")


            minimizedPdbBlock = self.igor.pose2str()

            if outdir:
                if molId is None:
                    molId = mol.GetProp("_Name")
                assert molId is not None,  "Error, molId required to save minimized structure"
                holo_fname = os.path.join(outdir, "{molId}.holo_minimised.pdb".format(molId=molId))
                with open(holo_fname, "w") as f:
                    f.write( minimizedPdbBlock )

                holo_fname = os.path.join(outdir, "{molId}.holo_unminimised.pdb".format(molId=molId))
                with open(holo_fname, "w") as f:
                    f.write( unminimizedPdbBlock )

                lig_fname = os.path.join(outdir, "{molId}.minimised.mol".format(molId=molId))
                Chem.MolToMolFile(ligand, lig_fname)

                md_fname = os.path.join(outdir, "{molId}.metadata.json".format(molId=molId))
                with open(md_fname, "w") as f:
                    json.dump(metadata, f)

            return ligand, metadata




def test():

    mol = Chem.Mol(b"\xef\xbe\xad\xde\x00\x00\x00\x00\x0c\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x0b\x00\x00\x00\x0b\x00\x00\x00\x80\x01\x07\x00h\x00\x00\x00\x03\x01\x02\x06\x00(\x00\x00\x00\x03\x04\x08\x00(\x00\x00\x00\x03\x02\x06\x00`\x00\x00\x00\x02\x02\x06\x00`\x00\x00\x00\x02\x02\x06@(\x00\x00\x00\x03\x04\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x06@h\x00\x00\x00\x03\x03\x01\x0b\x01\x00 \x02\x01(\x02\x03\x01\x00\x04\x03\x00\x05\x04\x00\x06\x05h\x0c\x07\x06h\x0c\x08\x07h\x0c\t\x08h\x0c\n\th\x0c\n\x05h\x0c\x17\x01\x00\x00\x00\x01\x00\x00\x00\x00\x0b1\x08\n\xc1\x98n*@\x8dW\x9a\xc2B`\x07\xc1^\xba\xa9?\xa40\x9a\xc2\xc5 \x16\xc1D\x8b\x0c?\x89\xc1\x99\xc2\xc3\xf5\xe0\xc0\x08\xac\\?\xa4p\x9a\xc2\xd3M\xce\xc0\xfe\xd4x?q=\x9d\xc2h\x91\xed\xc0\x8d\x97\x9e?^z\x9f\xc2-\xb2\x05\xc1\x85\xeb\x91>\x83@\xa0\xc2\x1f\x85\x13\xc1\xdd$\x06?\xb2]\xa2\xc2'1\x12\xc1\x1f\x85\xdb?\x06\xc1\xa3\xc2)\\\x03\xc1\\\x8f*@\xa2\x05\xa3\xc2=\n\xeb\xc0\x0c\x02\x1b@\xf6\xe8\xa0\xc2\x16")
    templateFname = "/home/ruben/oxford/myProjects/diamondCovid/data/nsp13/aligned/nsp13-x0176_0B/nsp13-x0176_0B_apo-desolv.pdb"
    wdir = "./output"
    results = MinimizePDBComplex_pyrosetta(templateFname).minimize( mol, wdir)
    print(results)

if __name__ == "__main__":

    print("trying minimize")
    test()

    '''

python -m fragmenstein.protocols.steps.minimizePdbComplex_pyrosetta

    '''