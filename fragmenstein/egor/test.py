# THIS IS BASICALLY A DEMO USAGE.

from . import Egor, pyrosetta

def test():
    acl = Egor.from_pdbfile(pdbfile='output/PAU-WEI-b9b-8_NIT2/pre_PAU-WEI-b9b-8_NIT2.pdb',
                            params_file='output/PAU-WEI-b9b-8_NIT2/PAU-WEI-b9b-8_NIT2.params',
                            constraint_file='cysBound.noCY.cst')
    acl.coordinate_constraint = 100
    print('initial')
    print(acl.ligand_score())
    print(acl.ligand_residue)  # AtomPair SG 145A HE2 41A HARMONIC 1.8 0.2
    print(acl.key_residues)
    r = acl.pose.pdb_info().pdb2pose(res=41, chain='A')
    MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
    MutateResidue(target=r, new_res='HIS_D').apply(acl.pose)
    MutateResidue(target=r, new_res='HIS').apply(acl.pose)
    pymover = pyrosetta.PyMOLMover()
    pymover.pymol_name(f'initial')
    native = acl.pose.clone()
    pymover.apply(native)
    acl.pose = native.clone()
    acl.minimise(10)
    pymover.pymol_name(f'relaxed')
    pymover.apply(acl.pose)
    print('cartesian relaxed')
    print(acl.ligand_score())
    acl.pose.dump_pdb('egor_test.pdb')