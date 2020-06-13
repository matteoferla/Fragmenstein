from ._victor_base_mixin import _VictorBaseMixin
from ..core import Fragmenstein
from ..igor import Igor
from ..m_rmsd import mRSMD
from ..rectifier import Rectifier
from typing import List, Optional, Dict, Union, Callable
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
import time

class _VictorAutomergeMixin(_VictorBaseMixin):

    @classmethod
    def combine(cls,
                hits: List[Chem.Mol],
                pdb_filename: str,
                ligand_resn: str = 'LIG',
                ligand_resi: Union[int, str] = '1B',
                covalent_resn: str = 'CYS',  # no other option is accepted.
                covalent_resi: Optional[Union[int, str]] = None,
                extra_constraint: Union[str] = None,
                pose_fx: Optional[Callable] = None,
                atomnames: Optional[Dict[int, str]] = None
                ):
        self = cls.__new__(cls)
        self.long_name = '-'.join([h.GetProp('_Name') for h in hits])
        self.apo_pdbblock = open(pdb_filename).read()
        self.hits = hits
        self.ligand_resn = ligand_resn.upper()
        self.ligand_resi = ligand_resi
        self.covalent_resn = covalent_resn.upper()
        self.covalent_resi = covalent_resi
        self.atomnames = atomnames
        self.extra_constraint = extra_constraint
        self.pose_fx = pose_fx
        # these are calculated
        self.is_covalent = None
        self.params = None
        self.mol = None
        self.constraint = None
        self.fragmenstein = None
        self.unminimised_pdbblock = None
        self.igor = None
        self.minimised_pdbblock = None
        self.minimised_mol = None
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.tick = time.time()
        self.tock = float('inf')
        self._safely_do(execute=self._combine_main,
                        resolve=self._resolve,
                        reject=self._reject)
        return self

    def _combine_main(self):
        self.fragmenstein = Fragmenstein(mol=Chem.Mol(),
                hits=[],
                attachment=None,
                merging_mode='off')
        # collapse hits
        self.fragmenstein.hits = [self.fragmenstein.collapse_ring(h) for h in self.hits]
        # merge!
        self.fragmenstein.scaffold = self.fragmenstein.merge_hits()
        self.journal.debug(f'{self.long_name} - Merged')
        self.fragmenstein.positioned_mol = self.fragmenstein.expand_ring(self.fragmenstein.scaffold, bonded_as_original=False)
        self.journal.debug(f'{self.long_name} - Expanded')
        self.fragmenstein.positioned_mol = Rectifier(self.fragmenstein.positioned_mol).mol
        # the origins are obscured because of the collapsing...
        self.fragmenstein.guess_origins(self.fragmenstein.positioned_mol, self.hits)
        self.fragmenstein.positioned_mol.SetProp('_Name', self.long_name)
        self.mol = self.fragmenstein.positioned_mol
        self.journal.debug(f'{self.long_name} - Rectified')
        self.smiles = Chem.MolToSmiles(self.mol)
        if self.fragmenstein_debug_draw:
            picture = Chem.CombineMols(Chem.CombineMols(self.hits[0], self.hits[1]), self.fragmenstein.positioned_mol)
            AllChem.Compute2DCoords(picture)
            self.fragmenstein.draw_nicely(picture)
        # making folder.
        self._make_output_folder()
        # paramterise
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        self.params = Params.load_mol(self.mol, name=self.ligand_resn)
        self.params.NAME = self.ligand_resn # force it.
        self.params.fix_mol()
        self.params.polish_mol()
        self.mol = Chem.Mol(self.params.mol)
        self.fragmenstein.positioned_mol = Chem.Mol(self.mol)
        self.params.mol = AllChem.AddHs(self.params.mol)
        AllChem.EmbedMolecule(self.params.mol)
        AllChem.MMFFOptimizeMolecule(self.params.mol)
        AllChem.ComputeGasteigerCharges(self.params.mol)
        # those Hs lack correct names
        self.params.fix_mol()
        self.params.convert_mol()
        self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        self.params.CHI.data = []  # TODO fix chi
        self._log_warnings()
        # get constraint
        # todo enable covalent
        ###self.constraint = self._get_constraint(self.extra_constraint)
        ####attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self.constraint = Constraints.mock()
        self.constraint.custom_constraint += self._make_coordinate_constraints_for_unnovels()
        self._log_warnings()
        self.post_params_step()
        self.fragmenstein_merging_mode = 'full'
        self.unminimised_pdbblock = self._place_fragmenstein()
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.unbound_pose = self._check_params()
        self._checkpoint_alpha()
        self._checkpoint_bravo()
        self.igor = Igor.from_pdbblock(pdbblock=self.unminimised_pdbblock,
                                       params_file=params_file,
                                       constraint_file=constraint_file,
                                       ligand_residue=self.ligand_resi,
                                       key_residues=[self.covalent_resi])
        # user custom code.
        if self.pose_fx is not None:
            self.journal.debug(f'{self.long_name} - running custom pose mod.')
            self.pose_fx(self.igor.pose)
        else:
            self.pose_mod_step()
        # storing a roundtrip
        self.unminimised_pdbblock = self.igor.pose2str()
        # minimise until the ddG is negative.
        ddG = self.reanimate()
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {ddG} kcal/mol {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        self.journal.debug(f'{self.long_name} - Completed')

    def _make_coordinate_constraints_for_unnovels(self):
        lines = []
        conf = self.fragmenstein.positioned_mol.GetConformer()
        for i, atom in enumerate(self.fragmenstein.positioned_mol.GetAtoms()):
            if atom.GetSymbol() == '*':
                continue
            elif atom.HasProp('_Novel') and atom.GetBoolProp('_Novel'):
                continue # novels
            else:
                pos = conf.GetAtomPosition(i)
                fxn = f'HARMONIC 0 1' # the other do not make sense here.
                lines.append(f'CoordinateConstraint {atom.GetPDBResidueInfo().GetName()} {self.ligand_resi} ' + \
                             f'CA {self.covalent_resi} ' + \
                             f'{pos.x} {pos.y} {pos.z} {fxn}\n')
        return ''.join(lines)