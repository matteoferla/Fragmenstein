import re, string, functools
from typing import Optional, Set

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_to_params import Params, Constraints
from ._victor_journal import _VictorJournal
from .minimalPDB import MinimalPDBParser
from ..monster._ff import MinizationOutcome


class _VictorPlonk(_VictorJournal):
    # the following is here as opposed to Monster, because it requires the template, ligand connections etc.
    # the stupid name "plonk" is to distinguish it from place and placement, which have a different meaning in Monster.

    def _get_LINK_record(self):
        if self.is_covalent:
            # get correct chain names.
            l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups()
            if self.covalent_resi:
                p_resi, p_chain = re.match('(\d+)(\D?)', str(self.covalent_resi)).groups()
            else:
                p_resi, p_chain = None, None
            if not p_chain:
                p_chain = 'A'
            if not l_chain:
                l_chain = 'B'
            # get the cx atom name
            assert len(self.params.CONNECT) != 0, 'No connection atom found in params file.'
            cx = self.params.pad_name(self.params.CONNECT[0].atom_name)
            # TODO the SG connection is hardcoded.
            return f'LINK         SG  {self.covalent_resn} {p_chain} {p_resi: >3}                ' + \
                f'{cx} {self.ligand_resn} {l_chain} {l_resi: >3}     1555   1555  1.8\n'
        else:
            return ''

    def _correct_ligand_info(self, mol: Optional[Chem.Mol] = None) -> Chem.Mol:
        """
        Corrects in place the given mol based on self.ligand_resi.
        If none provided it assumes ``self.mol`` (which in most cases is ``self.monster.positioned_mol``)
        Correcting the serial unfortunately does not do anything.

        Called by ``._plonk_monster_in_structure_minimal`` (a route of ``._plonk_monster_in_structure``)
        """
        i = 0
        if mol is None:
            mol = self.mol
        l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups()  # TODO improve ligand_resi
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is None:
                i += 1
                self.journal.warning(f'The atom #{atom.GetIdx()} has no PDB information in RDKit')
                info = Chem.AtomPDBResidueInfo(atomName=f'X{i: <3}')
            info.SetResidueNumber(int(l_resi))
            info.SetChainId(l_chain)
            info.SetIsHeteroAtom(True)
            info.SetOccupancy(1.)
            info.SetResidueName(self.ligand_resn)
        return mol

    def _plonk_monster_in_structure(self, prepped_mol: Optional[Chem.Mol]=None, use_pymol=False):
        """
        Plonks the molecule in the structure. This is most likely self.mol
        ``_correct_ligand_info`` will accept Chem.Mol or None (it will use self.mol).
        As the preminimisation is all that happens in Wictor, the prepped_mol is passed.
        """
        # this is in place on ``prepped_mol`` or ``self.mol`` (None case).
        # Passing its returned value to prepped_mol will result in minimisation being skipped (as its prepped already)
        self._correct_ligand_info(prepped_mol)
        self.journal.debug(f'{self.long_name} - placing monster in structure')
        if use_pymol:
            return self._plonk_monster_in_structure_pymol()
        else:
            # _plonk_monster_in_structure_raw does no corrections.
            return self._plonk_monster_in_structure_minimal(prepped_mol)

    @functools.cached_property
    def preminimized_mol(self) -> Chem.Mol:
        """
        This cached property is the preminimised molecule.
        It extracts the neighbourhood (``monster.get_neighborhood``) and
        minimises the molecule (``monster.mmff_minimize``)

        This method is called by the plonking into structure methods.
        Not "positioning" as intended by ``monster`` is done.
        """
        mol = Chem.Mol(self.mol)
        if self.monster_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising monster (MMFF)')
            if self.settings.get('ff_use_neighborhood', True):
                neighborhood = self.monster.get_neighborhood(self.apo_pdbblock,
                                                             cutoff=self.settings['ff_neighborhood'],
                                                             addHs=True)
            else:
                neighborhood = None
            # ff_max_displacement = float('nan') for fixed mode
            min_result: MinizationOutcome = self.monster.mmff_minimize(mol,
                                                    neighborhood=neighborhood,
                                                    ff_max_displacement=float(self.settings.get('ff_max_displacement', 0.)),
                                                    ff_constraint=int(self.settings.get('ff_constraint', 10)),
                                                   allow_lax=bool(self.settings['ff_prevent_cis']),
                                                   prevent_cis=bool( self.settings['ff_prevent_cis'])
                                                   )
            return min_result.mol
        else:
            return mol


    @functools.cached_property
    def preminimized_undummied_mol(self) -> Chem.Mol:
        """
        See ``preminimized_mol``. This strips the dummy atoms from the preminimised molecule.
        """
        return AllChem.DeleteSubstructs(self.preminimized_mol, Chem.MolFromSmiles('*'))

    def _plonk_monster_in_structure_minimal(self, prepped_mol: Optional[Chem.Mol]=None) -> str:
        """
        Plonks the molecule in the structure without using pymol.
        Uses a custom miniparser. see minimalPDB.MinimalPDBParser

        :return:
        """
        # ----- load
        if prepped_mol is None:
            mol = self.preminimized_undummied_mol
        else:  # this was passed by the user, fixing dummies just in case
            mol = AllChem.DeleteSubstructs(prepped_mol, Chem.MolFromSmiles('*'))
        pdbdata = MinimalPDBParser(self.apo_pdbblock, remove_other_hetatms=self.remove_other_hetatms,
                                   ligname=self.ligand_resn)
        moldata = MinimalPDBParser(Chem.MolToPDBBlock(mol))
        # ------- covalent fix
        if self.is_covalent:
            pdbdata.headers.append(self._get_LINK_record())
        # ------- assertions
        l_resi, l_chain = re.match('(\d+)(\D?)', str(self.ligand_resi)).groups()  # TODO improve ligand_resi
        if pdbdata.has_residue_index(index=int(l_resi), chain=l_chain):
            raise ValueError(f'Residue {self.ligand_resi} already exists in structure')
        elif pdbdata.has_residue_name(self.ligand_resn):
            raise ValueError(f'Residue {self.ligand_resn} already exists in structure')
        # -------- append
        pdbdata.append(moldata)  # fixes offsets in ATOM/HETATM and CONECT lines.
        return str(pdbdata)

    def _correct_covalent_resi(self):
        """
        An unresolved issue is that covalent_resi acts both as a covalent residue and the reference residue.
        This corrects for the case there is no covalent_resi
        """
        pdbdata = MinimalPDBParser(self.apo_pdbblock,
                                   remove_other_hetatms=self.remove_other_hetatms,
                                   ligname=self.ligand_resn)
        entry = pdbdata.coordinates[0]
        if self.covalent_resi is None:
            self.covalent_resi = f'{pdbdata.get_residue_index(entry)}{pdbdata.get_chain(entry)}'
        else:
            p_resi, p_chain = re.match('(\d+)(\D?)', str(self.covalent_resi)).groups()
            if not pdbdata.has_residue_index(int(p_resi), p_chain):
                self.covalent_resi = f'{pdbdata.get_residue_index(entry)}{pdbdata.get_chain(entry)}'

    def _get_empty_resi(self) -> str:
        """
        return the first empty chain basically.
        """
        pdbdata = MinimalPDBParser(self.apo_pdbblock,
                                   remove_other_hetatms=False)
        chains: Set[str] = {pdbdata.get_chain(entry) for entry in pdbdata.coordinates}
        missing = sorted(set(string.ascii_uppercase).difference(chains))
        return f'1{missing[0]}'

    def _plonk_monster_in_structure_raw(self):
        """
        Plonks the molecule in the structure without using pymol.
        Not "positioning" as intended by ``monster`` is done.
        Opening a PDB in RDKit is doable but gets exponentially slow with chain length

        :return:
        """
        mol = self.preminimized_undummied_mol
        mol_block = Chem.MolToPDBBlock(mol)
        return '\n'.join([self._get_LINK_record().strip(),
                          self.apo_pdbblock.strip(),
                          mol_block
                          ]).strip()

    def _plonk_monster_in_structure_pymol(self):
        """
        Plonks the molecule in the structure using pymol.
        Not "positioning" as intended by ``monster`` is done.

        :return:
        """
        import pymol2
        mol = self.preminimized_undummied_mol
        with pymol2.PyMOL() as pymol:
            pymol.cmd.read_pdbstr(self.apo_pdbblock, 'apo')
            pos_mol = Chem.MolToPDBBlock(mol)
            # pymol.cmd.read_pdbstr(pos_mol, 'scaffold')
            pymol.cmd.remove('name R')  # no dummy atoms!
            for c in self._connected_names:
                pymol.cmd.remove(f'name {c}')  # no conns
            pymol.cmd.remove('resn UNL')  # no unmatched stuff.
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return self._get_LINK_record() + pdbblock
