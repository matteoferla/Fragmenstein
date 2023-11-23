from ._utility import _MonsterUtil
from rdkit import Chem
from rdkit.Chem import AllChem, rdqueries, rdMolAlign
from typing import Optional, List, Union, Tuple, Dict
from warnings import warn
from dataclasses import dataclass
import numpy as np
import numpy.typing as npt


@dataclass
class MinizationOutcome:
    success: bool
    mol: Chem.Mol
    ideal: Chem.Mol
    U_pre: float = float('nan')
    U_post: float = float('nan')
    delta: float = float('nan')


class _MonsterFF(_MonsterUtil):

    def mmff_minimize(self,
                      mol: Optional[Chem.Mol] = None,
                      neighborhood: Union[Chem.Mol, None] = None,
                      ff_max_displacement: float = 0.,
                      ff_constraint: int = 10,
                      ff_max_iterations: int=200,
                      ff_cutoff: float = 100.,
                      allow_lax: bool = True) -> MinizationOutcome:
        """
        Minimises a mol, or self.positioned_mol if not provided, with MMFF constrained to ff_max_displacement Å.
        Gets called by Victor if the flag .monster_mmff_minimisation is true during PDB template construction.

        :param mol: Molecule to minimise. If None, self.positioned_mol is used.
        :param neighborhood: Protein neighboorhood (ignored if None)
        :param ff_max_displacement: Distance threshold (Å) for atomic positions mapped to hits for  MMFF constrains.
                            if NaN then fixed point constraints (no movement) are used.
                            This is passed as maxDispl to MMFFAddPositionConstraint.
        :param ff_constraint: Force constant for MMFF constraints.
        :param ff_cutoff: kcal/mol diff value to consider a failed minimisation.
        :param allow_lax: If True and the minimisation fails, the constraints are halved and the minimisation is rerun.
        :return: None

        Note that most methods calling this via Victor
        now use its ``.settings.py['ff_max_displacement']`` and ``.settings.py['ff_constraint']``
        and do not use the defaults.
        """
        # ## input fixes
        if mol is None and self.positioned_mol is None:
            raise ValueError('No valid molecule')
        elif mol is None:
            mol = self.positioned_mol
        else:
            pass  # mol is fine
        # ## prep
        success: bool
        fixed_mode = str(ff_max_displacement).lower() == 'nan'
        mol = Chem.Mol(mol)
        # protect
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_IsDummy')):
            atom.SetAtomicNum(8)
        combo, fixed_idxs = self._prep_combined(mol, neighborhood)
        ideal: Chem.Mol = self.make_ideal_mol(mol, ff_minimise=True)
        ideal_E: float = ideal.GetDoubleProp('Energy')
        # ## Start FF
        p = AllChem.MMFFGetMoleculeProperties(combo, 'MMFF94')
        if p is None:
            self.journal.error(f'MMFF cannot work on a molecule that has errors!')
            return MinizationOutcome(success=False, mol=mol, ideal=ideal)
        ff = AllChem.MMFFGetMoleculeForceField(combo, p, ignoreInterfragInteractions=False)
        if ff is None:
            return MinizationOutcome(success=False, mol=mol, ideal=ideal,
                                     U_post=float('nan'), U_pre=float('nan'), delta=0)
        # restrain
        restrained = []
        # mol not combo here:
        conserved: List[Chem.Atom] = list(mol.GetAtomsMatchingQuery(rdqueries.HasPropQueryAtom('_Novel', negate=True)))
        # weird corner case
        if len(conserved) == mol.GetNumAtoms() and fixed_mode:
            ff.Initialize()
            dU: float = ff.CalcEnergy()
            self.journal.warning('No novel atoms found in fixed_mode (ff_max_displacement == NaN), ' + \
                                 'this is probably a mistake')
            # nothing to do...
            return MinizationOutcome(success=True, mol=mol, ideal=ideal, U_post=dU, U_pre=dU, delta=0)
        # constrain or freeze
        for atom in conserved:
            i = atom.GetIdx()
            if atom.GetAtomicNum() == 1:
                # let hydrogens move
                continue
            elif fixed_mode:
                ff.AddFixedPoint(i)
            else:
                # https://github.com/rdkit/rdkit/blob/115317f43e3bdfd73673ca0e4c6b4035aa26a034/Code/ForceField/UFF/PositionConstraint.cpp#L35
                ff.MMFFAddPositionConstraint(i, maxDispl=ff_max_displacement, forceConstant=ff_constraint)
            restrained.append(i)
        # constrain dummy atoms
        for atom in mol.GetAtomsMatchingQuery(rdqueries.HasPropQueryAtom('_IsDummy')):
            i = atom.GetIdx()
            ff.MMFFAddPositionConstraint(i, maxDispl=0, forceConstant=ff_constraint * 5)
            restrained.append(i)
        for i in fixed_idxs:  # neighborhood is frozen
            ff.AddFixedPoint(i)
        # ## Minimize
        try:
            dG_pre = ff.CalcEnergy()
            dG_post = dG_pre
            previous_dG = float('inf')
            # this is a bit of a hack, but it works to make sure its not a flipped plateau-like local minima
            while previous_dG - dG_post > 0.5:
                previous_dG = dG_post
                m = ff.Minimize(maxIts=ff_max_iterations)
                dG_post = ff.CalcEnergy()
                if m == -1:
                    break
            if m == -1:
                self.journal.error('MMFF Minisation could not be started')
                success = False
            elif m == 0:
                self.journal.info('MMFF Minisation was successful')
                success = True
            elif m == 1:
                self.journal.info('MMFF Minisation was run, but the minimisation was not unsuccessful')
                success = False
            else:
                self.journal.critical("Iä! Iä! Cthulhu fhtagn! Ph'nglui mglw'nafh Cthulhu R'lyeh wgah'nagl fhtagn")
                success = False
        except RuntimeError as error:
            self.journal.info(f'MMFF minimisation failed {error.__class__.__name__}: {error}')
            return MinizationOutcome(success=False, mol=mol, ideal=ideal)
        # extract
        new_mol = self.extract_from_neighborhood(combo)
        ligand_E: float = self.MMFF_score(new_mol, delta=False)
        new_mol.SetDoubleProp('Energy', ligand_E)
        # check
        if ligand_E - ideal_E > abs(ff_cutoff):
            success = False  # damn
        if not success and allow_lax:
            self.journal.debug(f'MMFF minimisation failed, trying again with lax constraint {ff_constraint // 5}')
            return self.mmff_minimize(mol,
                                      neighborhood=neighborhood,
                                      ff_max_displacement=ff_max_displacement,
                                      ff_constraint=ff_constraint // 5,
                                      allow_lax=False)
        # deprotect
        for atom in new_mol.GetAtomsMatchingQuery(Chem.rdqueries.HasPropQueryAtom('_IsDummy')):
            atom.SetAtomicNum(0)
        # prevent drift:
        #rdMolAlign.AlignMol(new_mol, mol, atomMap=list(zip(restrained, restrained)))
        self.journal.info(f'MMFF minimisation: {dG_pre:.2f} -> {dG_post:.2f} kcal/mol '\
                           f'w/ {rdMolAlign.CalcRMS(new_mol, mol)}Å RMSD at '\
                           f'max displacement={ff_max_displacement} & constraint={ff_constraint}'
                           )
        return MinizationOutcome(success=success,
                                 mol=new_mol,
                                 ideal=ideal,
                                 U_post=dG_post,
                                 U_pre=dG_pre,
                                 delta=dG_post - dG_pre)

    def _prep_combined(self, mol, neighborhood) -> Tuple[Chem.Mol, List[int]]:
        # ## protect (DummyMasker could be used here)
        for atom in mol.GetAtomsMatchingQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(0)):
            atom.SetBoolProp('_IsDummy', True)
            atom.SetAtomicNum(16)
        Chem.SanitizeMol(mol)
        # ## Combine with neighborhood
        if neighborhood is not None:
            Chem.SanitizeMol(neighborhood)
            hydroneighborhood = AllChem.AddHs(neighborhood, addCoords=True)
            combo: Chem.Mol = Chem.CombineMols(mol, hydroneighborhood)
            Chem.SanitizeMol(combo)
            fixed_idxs: List[int] = list(range(mol.GetNumAtoms(), combo.GetNumAtoms()))
        else:
            combo = Chem.Mol(mol)
            fixed_idxs: List[int] = []
        self.journal.debug(f'Combined molecule has {combo.GetNumAtoms()} atoms, {fixed_idxs} fixed')
        return combo, fixed_idxs

    def MMFF_score(self, mol: Optional[Chem.Mol] = None, delta: bool = False, mode: str = 'MMFF') -> float:
        """
        Merck force field. Chosen over Universal for no reason at all.

        :param mol: ligand
        :type mol: ``Chem.Mol`` optional. If absent extracts from pose.
        :param delta: report difference from unbound (minimized)
        :type delta: bool
        :param mode: 'MMFF' or 'UFF'
        :type mode: str
        :return: kcal/mol
        :rtype: float

        :warning: This was moved out of Igor. Victor has the method for calling it with igor.mol_from_pose
        """
        if mol is None:
            mol = self.positioned_mol
        try:
            mol = Chem.Mol(mol)  # copy!
            if mode == 'UFF':
                ff = AllChem.UFFGetMoleculeForceField(mol)
            elif mode == 'MMFF':
                p = AllChem.MMFFGetMoleculeProperties(mol, 'MMFF94')
                ff = AllChem.MMFFGetMoleculeForceField(mol, p)
            else:
                raise ValueError(f'Unknown mode: {mode} (choice: MMFF or UFF)')
            ff.Initialize()
            # print(f'MMFF: {ff.CalcEnergy()} kcal/mol')
            if delta:
                pre = ff.CalcEnergy()
                ff.Minimize()
                post = ff.CalcEnergy()
                return pre - post
            else:
                return ff.CalcEnergy()
        except RuntimeError as err:
            self.journal.warning(f'{err.__class__.__name__}: {err} (It is generally due to bad sanitisation)')
            return float('nan')

    @classmethod
    def get_close_indices(cls, query: Chem.Mol, target: Chem.Mol, cutoff: float = 5.) -> List[int]:
        """
        Give an rdkit Chem.Mol ``query`` get the atom idices of ``target`` that are with ``cutoff`` Å.
        """
        combo = Chem.CombineMols(target, query)
        distances: npt.NDArray[np.float64] = AllChem.Get3DDistanceMatrix(combo)
        query2target_dist: npt.NDArray[np.float64] = distances[
            slice(target.GetNumAtoms(), combo.GetNumAtoms(), 1), slice(0, target.GetNumAtoms(), 1)].min(axis=0)
        neighbors: npt.NDArray[np.int64] = np.where(query2target_dist <= cutoff)[0]
        return list(map(int, neighbors))

    @classmethod
    def _get_aromatic_neighbors(cls, atom, accounted):
        neighbor: Chem.Atom
        for neighbor in atom.GetNeighbors():
            neigh_idx: int = neighbor.GetIdx()
            if neighbor.GetIsAromatic() and neigh_idx not in accounted:
                accounted.append(neigh_idx)
                return cls._get_aromatic_neighbors(neighbor, accounted)
        return accounted

    @classmethod
    def extract_atoms(cls, protein: Chem.Mol, keepers: List[int], expand_aromatics: bool = True) -> Chem.Mol:
        """
        Extract the given atom indices (``keepers``) from ``protein``.
        Expanding to full aromatic ring and copying conformers
        """
        pasteboard = Chem.RWMol()
        # ## Expand aromatic rings
        exkeepers = list(keepers)
        if expand_aromatics:
            for idx in keepers:
                atom: Chem.Atom = protein.GetAtomWithIdx(idx)
                if not atom.GetIsAromatic():
                    continue
                exkeepers = cls._get_aromatic_neighbors(atom, exkeepers)
        # ## Add atoms
        prot2paste: Dict[int, int] = {}
        for idx in exkeepers:
            atom: Chem.Atom = protein.GetAtomWithIdx(int(idx))  # no to np.int64
            if not expand_aromatics:
                atom.SetIsAromatic(False)
            prot2paste[idx] = pasteboard.AddAtom(atom)
        # ## Add bonds
        for prot_idx, paste_idx in prot2paste.items():
            atom: Chem.Atom = protein.GetAtomWithIdx(prot_idx)
            for prot_neighbor_idx in [n.GetIdx() for n in atom.GetNeighbors()]:
                if prot_neighbor_idx in prot2paste and prot_neighbor_idx > prot_idx:
                    prot_bond: Chem.Bond = protein.GetBondBetweenAtoms(prot_idx, prot_neighbor_idx)
                    paste_neighneighbor_idx = prot2paste[prot_neighbor_idx]
                    pasteboard.AddBond(paste_idx, paste_neighneighbor_idx,
                                       prot_bond.GetBondType() if expand_aromatics else Chem.BondType.SINGLE)

        pasteboard_conf = Chem.Conformer(len(keepers))
        positions: npt.NDArray = protein.GetConformer().GetPositions()
        for prot_idx, paste_idx in prot2paste.items():
            pasteboard_conf.SetAtomPosition(paste_idx, positions[prot_idx, :])
        pasteboard.AddConformer(pasteboard_conf)
        return pasteboard.GetMol()

    def get_neighborhood(self, apo_block: str, cutoff: float, mol: Optional[Chem.Mol] = None) -> Chem.Mol:
        if mol is None:
            mol = self.positioned_mol
        protein: Chem.Mol = Chem.MolFromPDBBlock(apo_block)
        neighbor_idxs: List[int] = self.get_close_indices(mol, protein, cutoff)
        neighborhood: Chem.Mol = self.extract_atoms(protein, neighbor_idxs)
        self.journal.debug(f'{cutoff}Å Neighborhood has {neighborhood.GetNumAtoms()} atoms')
        for atom in neighborhood.GetAtoms():
            atom.SetBoolProp('IsNeighborhood', True)
        return neighborhood

    def make_ideal_mol(self, mol: Optional[Chem.Mol]=None, ff_minimise: bool=False) -> Chem.Mol:
        if mol is None:
            mol = self.positioned_mol
        ideal = Chem.Mol(mol)
        ideal.SetDoubleProp('Energy', float('nan'))
        AllChem.EmbedMolecule(ideal)
        p = AllChem.MMFFGetMoleculeProperties(ideal, 'MMFF94')
        ff = AllChem.MMFFGetMoleculeForceField(ideal, p)
        ff.Initialize()
        if ff_minimise:
            ff.Minimize()
        energy = ff.CalcEnergy()
        ideal.SetDoubleProp('Energy', energy)
        return ideal

    def extract_from_neighborhood(self, mol: Chem.Mol) -> Chem.Mol:
        rwmol = Chem.RWMol(mol)
        rwmol.BeginBatchEdit()
        for atom in rwmol.GetAtoms():
            if atom.HasProp('IsNeighborhood'):
                rwmol.RemoveAtom(atom.GetIdx())
        rwmol.CommitBatchEdit()
        return Chem.GetMolFrags(rwmol.GetMol(), asMols=True)[0]