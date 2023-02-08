"""
Do note that monster contains a few visualisation functions.

"""

from io import StringIO
from fragmenstein.branding import divergent_colors
from IPython.display import display, HTML
from rdkit import Chem
import os, re
from typing import (Optional)
from ._victor_common import _VictorCommon
from ..display import ComponentViewer, MolNGLWidget

class _VictorShow(_VictorCommon):
    # partial move out of utility module

    def to_nglview(self, print_legend: bool = False) -> MolNGLWidget:
        """`
        generates a NGLWidget (``IPython.display.display`` will show it)
        with the compounds and the merged if present.
        To override the colours:
        The colours will be those in the ``Chem.Mol``'s property ``_color`` if present.

        Returns -> nv.NGLWidget
        """
        view, legend = self.monster._to_nglview_and_legend(show_positioned_mol=False)
        for molblock in (self.minimized_pdbblock, self.unminimized_pdbblock):
            if molblock is None:
                continue
            fh = StringIO(molblock)
            comp: ComponentViewer = view.add_component(fh, ext='pdb')  # noqa it's there.
            comp.add_representation('ball+stick', colorValue='white', multipleBond=True,
                                    sele=f'[{self.ligand_resn}]')
            # force it.
            comp.update_ball_and_stick(colorValue='white', multipleBond=True)
            legend += ' Fragmenstein monster (white)'
            view._js(f"""const comp = this.stage.compList[{len(self.hits)}]
                                 const target_sele = new NGL.Selection('[{self.ligand_resn}]');
                                 const radius = 5;
                                 const neigh_atoms = comp.structure.getAtomSetWithinSelection( target_sele, radius );
                                 const resi_atoms = comp.structure.getAtomSetWithinGroup( neigh_atoms );
                                 comp.addRepresentation( "line", {{sele: resi_atoms.toSeleString(), 
                                                                   colorValue: "gainsboro",
                                                                   multipleBond: true}} );
                                 comp.addRepresentation( "contact", {{sele: resi_atoms.toSeleString()}});
                            """)
            break
        else:
            pass
        if print_legend:
            display(HTML(legend))
        # async madness: disabled for now.
        # view.center(f'[{self.ligand_resn}]')
        return view


    # =================== save  ========================================================================================

    def make_pse(self,
                 filename: str = 'combo.pse',
                 extra_mols: Optional[Chem.Mol] = None):
        """
        Save a pse in the relevant folder. This is the Victor one.
        """
        assert '.pse' in filename, f'{filename} not .pse file'
        if extra_mols is None:
            extra_mols = []
        # ------------------------
        import pymol2
        with pymol2.PyMOL() as pymol:
            for hit in self.hits:
                hit_name = hit.GetProp('_Name')
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit), hit_name)
                if self.monster is None:
                    pymol.cmd.color('grey50', f'element C and {hit_name}')
                elif hit_name in self.monster.unmatched:
                    pymol.cmd.color('black', f'element C and {hit_name}')
                else:
                    pymol.cmd.color('white', f'element C and {hit_name}')
            if self.monster is not None and self.monster.positioned_mol is not None:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.monster.positioned_mol), 'placed')
                pymol.cmd.color('magenta', f'element C and placed')
                pymol.cmd.zoom('byres (placed expand 4)')
                pymol.cmd.show('line', 'byres (placed around 4)')
            if self.minimized_mol is not None:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(self.minimized_mol), 'minimised')
                pymol.cmd.color('green', f'element C and minimised')
            if self.minimized_pdbblock is not None:
                pymol.cmd.read_pdbstr(self.minimized_pdbblock, 'min_protein')
                pymol.cmd.color('gray50', f'element C and min_protein')
                pymol.cmd.hide('sticks', 'min_protein')
            if self.unminimized_pdbblock is not None:
                pymol.cmd.read_pdbstr(self.unminimized_pdbblock, 'unmin_protein')
                pymol.cmd.color('gray20', f'element C and unmin_protein')
                pymol.cmd.hide('sticks', 'unmin_protein')
                pymol.cmd.disable('unmin_protein')
            for mol in extra_mols:
                name = mol.GetProp('_Name')
                pymol.cmd.read_molstr(Chem.MolToMolBlock(mol, kekulize=False), name)
                pymol.cmd.color('magenta', f'{name} and name C*')
            pymol.cmd.save(os.path.join(self.work_path, self.long_name, filename))

    def make_steps_pse(self, filename: str = 'step.pse'):
        """
        Saves the steps in a pse file. For a more exhaustive file,
        use `make_pse` in Monster.
        """
        import pymol2
        assert '.pse' in filename, f'{filename} not .pse file'
        with pymol2.PyMOL() as pymol:
            for hit in self.hits:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(hit, kekulize=False), hit.GetProp('_Name'))
            for label, mod in self.modifications:
                pymol.cmd.read_molstr(Chem.MolToMolBlock(mod, kekulize=False), re.sub('[^\w_]', '_', label))
            pymol.cmd.save(os.path.join(self.work_path, self.long_name, filename))

    # ===== MF keeps calling these from victor...

    def draw_nicely(self, *args, **kwargs):
        self.journal.warning('The `draw_nicely` method is in monster. Please use Monster.draw_nicely')
        self.monster.draw_nicely(*args, **kwargs)

    def show_comparison(self, *args, **kwargs):
        self.journal.warning('The `show_comparison` method is in monster. Please use Monster.show_comparison')
        self.monster.show_comparison(*args, **kwargs)


