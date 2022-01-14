import json
import os

from rdkit import Chem

from ._victor_overridables import _VictorOverridables


class _VictorStore(_VictorOverridables):
    # _save_prerequisites is in VictorCommon

    def checkpoint(self):
        self._checkpoint_alpha()
        self._checkpoint_bravo()
        self._checkpoint_charlie()

    def _checkpoint_alpha(self):
        from ..igor.pyrosetta_import import pyrosetta
        self._log_warnings()
        # saving hits (without copying)
        for h, hit in enumerate(self.hits):
            if hit.HasProp("_Name") and hit.GetProp("_Name").strip():
                name = hit.GetProp("_Name")
            else:
                name = f'hit{h}'
            hfile = os.path.join(self.work_path, self.long_name, f'{name}.pdb')
            Chem.MolToPDBFile(hit, hfile)
            mfile = os.path.join(self.work_path, self.long_name, f'{name}.mol')
            Chem.MolToMolFile(hit, mfile, kekulize=False)
        # saving params template
        params_template_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_template.pdb')
        Chem.MolToPDBFile(self.params.mol, params_template_file)
        params_template_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_template.mol')
        Chem.MolToMolFile(self.params.mol, params_template_file)
        # checking all is in order
        ptest_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_test.pdb')
        self.unbound_pose.dump_pdb(ptest_file)
        pscore_file = os.path.join(self.work_path, self.long_name, self.long_name + '.params_test.score')
        scorefxn = pyrosetta.get_fa_scorefxn()
        with open(pscore_file, 'w') as w:
            w.write(str(scorefxn(self.unbound_pose)))
        self._log_warnings()

    def _checkpoint_bravo(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving mols from monster')
        # if self.monster.scaffold is not None:
        #     scaffold_file = os.path.join(self.work_path, self.long_name, self.long_name + '.scaffold.mol')
        #     Chem.MolToMolFile(self.monster.scaffold, scaffold_file, kekulize=False)
        #     if self.monster.scaffold.HasProp('parts'):
        #         disregard = json.loads(self.monster.scaffold.GetProp('parts'))
        #         self.journal.info(f'{self.long_name} - disregarded {disregard}')
        #     else:
        #         disregard = []
        # if self.monster.chimera is not None:
        #     chimera_file = os.path.join(self.work_path, self.long_name, self.long_name + '.chimera.mol')
        #     Chem.MolToMolFile(self.monster.chimera, chimera_file, kekulize=False)
        if self.monster.positioned_mol is not None:
            pos_file = os.path.join(self.work_path, self.long_name, self.long_name + '.positioned.mol')
            Chem.MolToMolFile(self.monster.positioned_mol, pos_file, kekulize=False)
        if self.monster.mol_options:
            opt_file = os.path.join(self.work_path, self.long_name, self.long_name + '.mol_options.sdf')
            writer = Chem.SDWriter(opt_file)
            writer.SetKekulize(False)
            for t in self.monster.mol_options:
                writer.write(t)
            writer.close()
        frag_file = os.path.join(self.work_path, self.long_name, self.long_name + '.monster.json')
        data = {'smiles': self.smiles,
                'origin': self.monster.origin_from_mol(self.monster.positioned_mol),
                'stdev': self.monster.stdev_from_mol(self.monster.positioned_mol)}
        # if disregard:
        #     data['disregard'] = disregard
        with open(frag_file, 'w') as w:
            json.dump(data, w)
        self._log_warnings()
        # unminimized_pdbblock will be saved by igor (round trip via pose)

    def _checkpoint_charlie(self):
        self._log_warnings()
        self.journal.debug(f'{self.long_name} - saving pose from igor')
        min_file = os.path.join(self.work_path, self.long_name, self.long_name + '.holo_minimised.pdb')
        self.igor.pose.dump_pdb(min_file)
        self.journal.debug(f'{self.long_name} - calculating Gibbs')
        # recover bonds
        lig_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.mol')
        Chem.MolToMolFile(self.minimized_mol, lig_file)
        score_file = os.path.join(self.work_path, self.long_name, self.long_name + '.minimised.json')
        with open(score_file, 'w') as w:
            json.dump({'Energy': self.energy_score,
                       'mRMSD': self.mrmsd.mrmsd,
                       'RMSDs': self.mrmsd.rmsds}, w)
        self._log_warnings()