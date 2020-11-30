# this method is kept here as it has high salvage potential...
# namely it predicts what the covalent points are etc.

class DeleteMe:
    @classmethod
    def validate(cls,
                  hits: List[Chem.Mol],
                  smiles: str,
                  followup_pdb_filename: str,
                  template_pdb_filename: str,
                  ligand_resn: str = 'LIG'):

        """
        Given a pdb file with the protein with ligand and a smiles string, score it against hits.
        Without changing the position of the ligand —bar for protonations.
        The purpose of this is to get a DDG and a RMSD of a given followup for validation purposes.
        It is a bit convoluted solely to have consistent atom names with the followup,
        as opposed to changing them afterwards —safer for any arbitrary test.
        The followup pdb is not used for fragmenstein placement as that would bias the test!

        :param hits:
        :param smiles:
        :param followup_pdb_filename: pdb file with followup
        :param template_pdb_filename: pdb file of apo structure to use for placement

        :param ligand_resn:
        :return:
        """
        raise ValueError('THIS MAY BE BIASED!')
        # ***** STORE *******
        # entry attributes
        self = cls.__new__(cls)
        self.adam_merging_mode = 'none_permissive'
        fname = long_name = os.path.split(followup_pdb_filename)[-1]
        long_name = os.path.splitext(fname)[0] + '_validation'
        self.long_name = self.slugify(long_name)
        self.pdb_filename = followup_pdb_filename # non standard!
        self.unminimised_pdbblock = None
        self.apo_pdbblock = open(template_pdb_filename).read()
        self.hits = hits
        self.ligand_resn = ligand_resn
        pdb = Chem.MolFromPDBFile(followup_pdb_filename)
        attachment, attachee = self.find_attachment(pdb, ligand_resn)
        if attachment is not None:
            self.is_covalent = True
            if '*' in smiles:
                self.smiles = smiles
            else:
                self.smiles = self.make_covalent(smiles)
        else:
            self.smiles = smiles
            self.is_covalent = False
            attachment, attachee = self.find_closest_to_ligand(pdb, ligand_resn)
        info = attachment.GetPDBResidueInfo()
        self.covalent_resn = info.GetResidueName()
        self.covalent_resi = str(info.GetResidueNumber() ) +info.GetChainId()
        info = attachee.GetPDBResidueInfo()
        self.ligand_resi = str(info.GetResidueNumber() ) +info.GetChainId()
        self.pose_fx = None
        # these are calculated
        self.params = None
        self.mol = None
        self.constraint = None
        self.extra_constraint = None
        self.adam = None
        self.igor = None
        self.minimised_pdbblock = None
        self.minimised_mol = None
        # this is unique to validate
        self.reference_mol = self.extract_mol(name='crystal',
                                              filepath=followup_pdb_filename,
                                              smiles=smiles,
                                              ligand_resn=ligand_resn)
        # buffers etc.
        self._warned = []
        self.energy_score = {'ligand_ref2015': {'total_score': float('nan')},
                             'unbound_ref2015': {'total_score': float('nan')}}
        self.mrmsd = mRSMD.mock()
        self.tick = time.time()
        self.tock = float('inf')
        # analyse
        self._safely_do(execute=self._vanalyse, resolve=self._resolve, reject=self._reject)
        return self

    def _vanalyse(self):
        # THIS IS A COPY PASTE EXCEPT FOR REANIMATE and Params!!
        # self._assert_inputs()
        # ***** PARAMS & CONSTRAINT *******
        self.journal.info(f'{self.long_name} - Starting work')
        self._log_warnings()
        # making folder.
        self._make_output_folder()
        # make params
        self.journal.debug(f'{self.long_name} - Starting parameterisation')
        self.params = Params.from_smiles_w_pdbfile(pdb_file=self.pdb_filename,
                                                   smiles=self.smiles, generic= False,
                                                   name=self.ligand_resn,
                                                   proximityBonding=False)
        self.journal.warning(f'{self.long_name} - CHI HAS BEEN DISABLED')
        self.params.CHI.data = []  # TODO fix chi
        self.mol = self.params.mol
        self._log_warnings()
        # get constraint
        self.constraint = self._get_constraint(self.extra_constraint)
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self._log_warnings()
        self.post_params_step()
        # ***** FRAGMENSTEIN *******
        # make adam
        attachment = self._get_attachment_from_pdbblock() if self.is_covalent else None
        self.journal.debug(f'{self.long_name} - Starting adam')
        self.adam = Adam(mol=self.params.mol,  # Chem.MolFromSmiles(self.smiles)
                                         hits=self.hits,
                                         attachment=attachment,
                                         merging_mode=self.adam_merging_mode,
                                         debug_draw=self.adam_debug_draw,
                                         average_position=self.adam_average_position)
        self.unminimised_pdbblock = self._place_adam()
        if self.adam_mmff_minisation:
            self.journal.debug(f'{self.long_name} - pre-minimising adam (MMFF)')
            self.adam.mmff_minimise(self.adam.positioned_mol)
        self.constraint.custom_constraint += self._make_coordinate_constraints()
        self._checkpoint_bravo()
        # save stuff
        params_file, holo_file, constraint_file = self._save_prerequisites()
        self.post_adam_step()
        self.unbound_pose = self.params.test()
        self._checkpoint_alpha()
        # ***** EGOR *******
        self.journal.debug(f'{self.long_name} - setting up Igor')
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
        # DO NOT DO ddG = self.reanimate()
        # ddG = self.quick_renamiate() # put igor to work!
        ddG = self.reanimate()
        self.minimised_pdbblock = self.igor.pose2str()
        self.post_igor_step()
        self.minimised_mol = self._fix_minimised()
        self.mrmsd = self._calculate_rmsd()
        self.journal.info(f'{self.long_name} - final score: {ddG} kcal/mol {self.mrmsd.mrmsd}.')
        self._checkpoint_charlie()
        # RMSD against self.reference_mol and docked_mol
        m = mRSMD.from_other_annotated_mols(self.minimised_mol, self.hits, self.adam.positioned_mol)

        docked_mol = self.dock()
        # RMSD again
        self.journal.debug(f'{self.long_name} - Completed')
