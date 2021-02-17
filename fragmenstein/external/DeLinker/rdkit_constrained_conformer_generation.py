#!/usr/bin/env python

import sys
import traceback
from typing import List, Union, Dict

import rdkit

try:
    from collections.abc import Iterable  # noqa
except ImportError:
    from collections import Iterable  # noqa


from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem
import numpy as np

from joblib import Parallel, delayed

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')

'''
Given a list of SMILES strings, generate 3D conformers in sdf format using RDKit.  
Energy minimizes and filters conformers to meet energy window and rms constraints.
Script modified from: https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py
In particular, constrained optimisation (with slack) added.
'''

#convert smiles to sdf
def getRMS(mol, c1,c2):
    (rms,trans) = AllChem.GetAlignmentTransform(mol,mol,c1,c2)
    return rms

def getRMS2(mol1, mol2):
    (rms,trans) = AllChem.GetAlignmentTransform(mol1,mol2)
    return rms


def RoughConstrainedEmbed(mol, core_with_exit, tdist=0.25, randomseed=2342,
                          getForceField=AllChem.UFFGetMoleculeForceField, **kwargs):
    """ generates an embedding of a molecule where part of the molecule
    is constrained to have particular coordinates
    Arguments
    :param mol: the molecule to embed
    :param core_with_exit: the molecule to use as a source of constraints. It must contain the exit vectors
    :param tdsit. Max distance constrain between mol and core_with_exit atoms
    :param randomSeed: (optional) seed for the random number generator
    """
    Chem.SanitizeMol(core_with_exit)

    # Renumber dummy atoms to end of core
    dummy_idx = []
    for atom in core_with_exit.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummy_idx.append(atom.GetIdx())
    sub_idx = list(range(core_with_exit.GetNumHeavyAtoms()+2))
    for idx in dummy_idx:
        sub_idx.remove(idx)
        sub_idx.append(idx)
    mol_range = list(range(core_with_exit.GetNumHeavyAtoms()+2))
    idx_to_add = list(set(mol_range).difference(set(sub_idx)))
    sub_idx.extend(idx_to_add)
    aligned_core_with_exit = Chem.rdmolops.RenumberAtoms(core_with_exit, sub_idx)

    # Match constrained substructure to full molecule
    du = Chem.MolFromSmiles('*')
    qp = Chem.AdjustQueryParameters()
    qp.makeDummiesQueries=True
    qlink = Chem.AdjustQueryProperties(aligned_core_with_exit,qp)
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    # from matplotlib import pyplot as plt; from rdkit.Chem import Draw
    # plt.imshow(Draw.MolsToGridImage([ qlink], molsPerRow=1));
    # plt.show()

    match = mol.GetSubstructMatch(qlink)

    if not match:
        print("mol: %s core_with_exit: %s qlink: %s"%(Chem.MolToSmiles(mol), Chem.MolToSmiles(core_with_exit), Chem.MolToSmiles(qlink)))
        raise ValueError("molecule doesn't match the core or exit vectors joined")

    # Generate conformer
    mol = Chem.AddHs(mol)
    ci = AllChem.EmbedMolecule(mol, randomSeed=randomseed, **kwargs)
    if ci < 0:
        raise ValueError('Could not embed molecule.')

    algMap = [(j, i) for i, j in enumerate(match) if i in list(range(aligned_core_with_exit.GetNumHeavyAtoms()))]
    # rotate the embedded conformation onto the core:
    aligned_core_without_exit = Chem.RemoveHs(AllChem.ReplaceSubstructs(aligned_core_with_exit,du,Chem.MolFromSmiles('[H]'),True)[0])
    rms = AllChem.AlignMol(mol, aligned_core_without_exit, atomMap=algMap)
    ff = getForceField(mol, confId=0)
    conf = aligned_core_without_exit.GetConformer()
    for i in range(aligned_core_without_exit.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
        ff.AddDistanceConstraint(pIdx, match[i], 0, tdist, 100.)
    ff.Initialize()
    n = 4
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
        n -= 1
    # realign
    rms = AllChem.AlignMol(mol, aligned_core_without_exit, atomMap=algMap)
    # get energy
    energy = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
    return mol, energy


def gen_confs_constrained(smiles_full_mols:List[List[rdkit.Chem.Mol]], sdffile: Union[str, None], reference_confs: List[List[rdkit.Chem.Mol]],
                          baseconfs=20, maxconfs=20, sample_mult=1, seed=42, rms_threshold=0.7, energy=10, 
                          verbose=False, mmff=False, nomin=False, etkdg=False, tdist=0.75,
                          smi_frags: rdkit.Chem.Mol=[], numcores=20, jpsettings=False, every=5)-> List[Dict[str, List[rdkit.Chem.Mol]]]:
    """ generates constrained conformations of molecules
    Arguments
    :param smiles_full_mols: list of lists of SMILES strings to generate conformers (one list containing all mols with same constrained substructure)
    :param sdffile: desired location of output file. If None, no output file will be written
    :param reference_confs: SD file containined conformations of constrained portion of molecule (with dummy/R atom type present for exit vectors)
    :param baseconfs: how many unconstrained conformations to generate (used for energy gap)
    :param maxconfs: maximum number of conformers to return
    :param sample_mult: multiplier for number of conformers to generate
    :param seed: random seed
    :param rms_threshold: minimum RMS between conformers
    :param energy: maximum energy gap to lowest energy conformer
    :param mmff: use MMFF forcefield
    :param nomin: True if do not minimise
    :param etkdg: True if use ETKDG
    :param tdist: Tether distance. Amount of slack in conformation of constained atoms (in Angstrom)
    :param smi_frags: list of SMILES strings of constrained portion of molecule (optional - only used for labelling purposes)
    :param numcores: number of CPU cores for conformer generation
    :param jpsettings: whether to use the sampling procedure described in Ebejer et al., 2012 (https://pubs.acs.org/doi/10.1021/ci2004658)
    :param every: how often to print out progress (set to a large number and set verbose to False to stop most outputs)
    :return results_conformers: a list of dicts(smile->conformers) that contains different conformers for the different input smiles
    """

    if sdffile:
      outf = open(sdffile,'w+')
      sdwriter = Chem.SDWriter(outf)
      if sdwriter is None:
          print("Could not open ".sdffile)
          sys.exit(-1)

    if etkdg and not AllChem.ETKDG:
        print("ETKDG does not appear to be implemented.  Please upgrade RDKit.")
        sys.exit(1)

    if smi_frags != []:
        if len(smiles_full_mols) != len(smi_frags):
            print("smiles_full_mols and smi_frags not equal in length")
            return None

    if jpsettings == True:
        rms_threshold = 0.35
        sample_mult = 1

    if not isinstance(reference_confs, Iterable):
        sdreader = Chem.SDMolSupplier(reference_confs)
    else:
        sdreader= reference_confs

    results_conformers=[]
    for count, (smis, reference) in enumerate(zip(smiles_full_mols, sdreader)):
        proposed_conformers = {}
        for smi in smis:
            proposed_conformers[smi]=[]
            name = smi
            if count % every == 0:
                print("\rProcessed: %d " % count, end='')

            pieces = smi.split('.')
            if len(pieces) > 1:
                smi = max(pieces, key=len) #take largest component by length
                print("Taking largest component: %s" % (smi))

            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                if verbose:
                    print("smi:", smi)
                try:
                    Chem.SanitizeMol(mol)
                    mol = Chem.AddHs(mol)
                    mol.SetProp("_Name", name)
                    if jpsettings == True:
                        rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
                        if rot_bonds <= 7:
                            maxconfs=50
                        elif rot_bonds >=8 and rot_bonds <= 12:
                            maxconfs=200
                        else:
                            maxconfs=300
                    if smi_frags != []:
                        mol.SetProp("_StartingPoint", smi_frags[count][0])

                    mol_no_confs = Chem.Mol(mol)

                    # Initially embed sample_mult*baseconfs times
                    if etkdg:
                        #cids = Chem.EmbedMultipleConfs(mol, numConfs=int(sample_mult*maxconfs), params=Chem.ETKDG())#, numThreads=0)
                        cids = AllChem.EmbedMultipleConfs(mol, numConfs=int(sample_mult*baseconfs), useExpTorsionAnglePrefs=True, useBasicKnowledge=True, randomSeed=seed, numThreads=numcores)
                    else:
                        cids = AllChem.EmbedMultipleConfs(mol, int(sample_mult*baseconfs),randomSeed=seed, numThreads=numcores)
                    if verbose:
                        print(len(cids),"conformers found")
                    # And minimise
                    cenergy = []            
                    if mmff:
                        converged_res = AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=numcores)
                        cenergy = [i[1] for i in converged_res]
                    elif not nomin and not mmff:
                        converged_res = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=numcores)
                        cenergy = [i[1] for i in converged_res]
                    else:       
                        for conf in cids:
                            #not passing confID only minimizes the first conformer
                            if nomin:
                                cenergy.append(conf)
                            elif mmff:
                                converged = AllChem.MMFFOptimizeMolecule(mol,confId=conf)
                                mp = AllChem.MMFFGetMoleculeProperties(mol)
                                cenergy.append(AllChem.MMFFGetMoleculeForceField(mol,mp,confId=conf).CalcEnergy()) 
                            else:
                                converged = not AllChem.UFFOptimizeMolecule(mol,confId=conf)
                                cenergy.append(AllChem.UFFGetMoleculeForceField(mol,confId=conf).CalcEnergy())
                            if verbose:
                                print("Convergence of conformer",conf,converged)
                
                    # Now do constrained embedding
                    backend= 'multiprocessing' #'multiprocessing'
                    with Parallel(n_jobs=numcores, backend=backend) as parallel:
                        np.random.seed(seed=seed)
                        seeds = [int(x) for x in np.random.randint(1,1000001, sample_mult*maxconfs)]
                        results = parallel(delayed(RoughConstrainedEmbed)(Chem.Mol(mol_no_confs), reference, tdist=tdist, randomseed=s) for s in seeds)
                        # results= list(map(lambda rseed: RoughConstrainedEmbed(Chem.Mol(mol_no_confs), reference, tdist=tdist, randomseed=rseed), seeds))
                        mols = [r[0] for r in results]
                        cenergy_con = [r[1] for r in results]

                    # Filter by RMS and energy
                    if len(cenergy)>0 and len(cenergy_con)>0:
                        mine = min(cenergy+cenergy_con)
                        # TEST
                        #mine = min(cenergy_con)

                    mol = Chem.RemoveHs(mol)
                    #sortedcids = sorted(cids,key = lambda cid: cenergy[cid])
                    if(rms_threshold == 0):
                        cnt = 0;
                        for conf_num, conf in enumerate(mols):
                            if(cnt >= maxconfs):
                                break
                            if(energy < 0) or cenergy_con[conf_num]-mine <= energy:
                                conf.SetProp("_Model", str(conf_num))
                                conf.SetProp("_StartingPoint", smi_frags[count][0])
                                conf.SetProp("_Name", smi)
                                if sdffile:
                                  sdwriter.write(conf)
                                cnt+=1
                                proposed_conformers[smi].append(conf)
                    else:
                        written = []
                        for conf_num, conf in enumerate(mols):
                            if len(written) >= maxconfs:
                                break
                            #check rmsd
                            passed = True
                            for seenconf in written:
                                rms = getRMS2(conf,seenconf) 
                                if(rms < rms_threshold) or (energy > 0 and cenergy_con[conf_num]-mine > energy):
                                    passed = False
                                    break
                            if(passed):
                                written.append(conf)
                                conf.SetProp("_Model", str(conf_num))
                                conf.SetProp("_StartingPoint", smi_frags[count][0])
                                conf.SetProp("_Name", smi)
                                if sdffile:
                                  sdwriter.write(conf)
                                proposed_conformers[smi].append(conf)
                except (KeyboardInterrupt, SystemExit):
                    raise                
                except Exception as e:
                    traceback.print_exc()
                    print("Exception",e)
            else:
                print("ERROR:",smi)

        results_conformers.append(proposed_conformers)

    if sdffile:
      sdwriter.close()
    return results_conformers

if __name__ == "__main__":
    exit()
