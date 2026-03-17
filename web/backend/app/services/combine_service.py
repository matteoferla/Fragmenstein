"""Wraps Laboratory.combine() for the web API."""

import logging
from pathlib import Path
from typing import List

from rdkit import Chem

from fragmenstein import Laboratory, Victor, Wictor
from fragmenstein.laboratory._base import binarize, unbinarize

from . import file_manager
from .result_serializer import save_dataframe

log = logging.getLogger(__name__)

VICTOR_TYPES = {
    "Victor": Victor,
    "Wictor": Wictor,
}

# Import optional variants
try:
    from fragmenstein import Quicktor
    VICTOR_TYPES["Quicktor"] = Quicktor
except ImportError:
    pass

try:
    from fragmenstein import OpenVictor
    VICTOR_TYPES["OpenVictor"] = OpenVictor
except ImportError:
    pass


class WebLaboratory(Laboratory):
    """Laboratory subclass that passes warhead_harmonisation and merging_mode
    through to Victor subprocess calls (the base class hardcodes defaults)."""

    _warhead_harmonisation: str = "first"
    _merging_mode: str = "expansion"

    def combine_subprocess(self, binary_hits: List[bytes]) -> dict:
        """Override to pass warhead_harmonisation to victor.combine()."""
        from fragmenstein.igor import pyrosetta

        if self.Victor.uses_pyrosetta:
            pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        tentative_name = 'UNKNOWN'
        try:
            self.journal.debug(f'Combining {len(binary_hits)} hits')
            hits: List[Chem.Mol] = [hit for hit in map(unbinarize, binary_hits) if hit is not None]
            assert len(hits) > 0, f'No valid hits ({len(binary_hits)} provided)'
            assert all([hit.GetNumAtoms() > 0 for hit in hits]), 'Some hits have no atoms!'
            if all([mol.HasProp('_Name') for mol in hits]):
                tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
                if tentative_name in self.blacklist:
                    raise ValueError(f'{tentative_name} is blacklisted')
            self.journal.debug(f'Using {self.Victor.__name__}')
            victor = self.Victor(hits=hits,
                                 pdb_block=self.pdbblock,
                                 ligand_resn='LIG',
                                 ligand_resi=self.ligand_resi,
                                 covalent_resi=self.covalent_resi,
                                 **self.settings)
            tentative_name = '-'.join([mol.GetProp('_Name') for mol in hits])
            if tentative_name in self.blacklist:
                raise ValueError(f'{tentative_name} is blacklisted')
            victor.monster_throw_on_discard = True
            victor.monster.throw_on_discard = True
            victor.combine(warhead_harmonisation=self._warhead_harmonisation)
            result: dict = victor.summarize()
            result['unmin_binary'] = binarize(victor.monster.positioned_mol)
            result['min_binary'] = binarize(victor.minimized_mol)
            result['hit_binaries'] = [binarize(h) for h in victor.hits]
            if self.run_plip:
                result.update(victor.get_plip_interactions())
            return result
        except KeyboardInterrupt as err:
            raise err
        except Exception as error:
            error_msg = f'{error.__class__.__name__} {error}'
            self.Victor.journal.critical(f'*** {error_msg} for {tentative_name}')
            return dict(error=error_msg, name=tentative_name)

    def place_subprocess(self, inputs) -> dict:
        """Override to pass merging_mode to victor.place()."""
        from fragmenstein.igor import pyrosetta

        name: str = inputs['name']
        smiles: str = inputs['smiles']
        if self.Victor.uses_pyrosetta:
            pyrosetta.distributed.maybe_init(extra_options=self.init_options)
        try:
            binary_hits = inputs['binary_hits']
            hits: List[Chem.Mol] = [hit for hit in map(unbinarize, binary_hits) if hit]
            assert len(hits) > 0, 'No valid hits!'
            self.journal.debug(f'Using {self.Victor.__name__}')
            victor = self.Victor(hits=hits,
                                 pdb_block=self.pdbblock,
                                 ligand_resn='LIG',
                                 ligand_resi=self.ligand_resi,
                                 covalent_resi=self.covalent_resi,
                                 **self.settings)
            victor.place(smiles, long_name=name, merging_mode=self._merging_mode)
            result: dict = {**dict(inputs), **victor.summarize()}
            result['unmin_binary'] = binarize(victor.monster.positioned_mol)
            result['min_binary'] = binarize(victor.minimized_mol)
            result['hit_binaries'] = [binarize(h) for h in victor.hits]
            if self.run_plip and victor.minimized_pdbblock:
                result.update(victor.get_plip_interactions())
            return result
        except KeyboardInterrupt as err:
            raise err
        except Exception as error:
            error_msg = f'{error.__class__.__name__} {error}'
            self.Victor.journal.critical(f'*** {error_msg} for {name}')
            return dict(error=error_msg, name=name)


def load_session_hits(session_id: str, hit_names: list[str] | None = None) -> list[Chem.Mol]:
    """Load all hit molecules for a session from the database."""
    from ..database import get_db

    with get_db() as db:
        if hit_names:
            placeholders = ",".join("?" * len(hit_names))
            rows = db.execute(
                f"SELECT name, mol_block FROM hits WHERE session_id = ? AND name IN ({placeholders})",
                [session_id] + hit_names,
            ).fetchall()
        else:
            rows = db.execute(
                "SELECT name, mol_block FROM hits WHERE session_id = ?",
                (session_id,),
            ).fetchall()

    mols = []
    for row in rows:
        mol = Chem.MolFromMolBlock(row["mol_block"], removeHs=False)
        if mol is not None:
            mol.SetProp("_Name", row["name"])
            mols.append(mol)
    return mols


def load_session_template(session_id: str) -> str:
    """Load the template PDB block for a session."""
    from ..database import get_db

    with get_db() as db:
        row = db.execute(
            "SELECT template_pdb FROM sessions WHERE id = ?",
            (session_id,),
        ).fetchone()

    if row is None or row["template_pdb"] is None:
        raise ValueError(f"No template uploaded for session {session_id}")
    return row["template_pdb"]


def run_combine(
    session_id: str,
    victor_type: str = "Wictor",
    n_cores: int = -1,
    timeout: int = 240,
    combination_size: int = 2,
    permute: bool = True,
    joining_cutoff: float = 5.0,
    quick_reanimation: bool = False,
    warhead_harmonisation: str = "first",
    run_plip: bool = False,
    covalent_resi: str | None = None,
    hit_names: list[str] | None = None,
) -> Path:
    """Run Laboratory.combine() and save results. Returns path to results file."""
    hits = load_session_hits(session_id, hit_names)
    if not hits:
        raise ValueError("No hits loaded for session")

    pdbblock = load_session_template(session_id)
    work_dir = file_manager.session_work_dir(session_id)

    # Resolve Victor type
    VictorClass = VICTOR_TYPES.get(victor_type, Wictor)

    # Configure Victor class attributes
    WebLaboratory.Victor = VictorClass
    WebLaboratory.Victor.work_path = str(work_dir)
    WebLaboratory.Victor.monster_throw_on_discard = True
    WebLaboratory.Victor.monster_joining_cutoff = joining_cutoff
    WebLaboratory.Victor.quick_reanimation = quick_reanimation
    WebLaboratory.Victor.error_to_catch = Exception

    lab = WebLaboratory(pdbblock=pdbblock, covalent_resi=covalent_resi, run_plip=run_plip)
    lab._warhead_harmonisation = warhead_harmonisation

    df = lab.combine(
        hits,
        n_cores=n_cores,
        timeout=timeout,
        combination_size=combination_size,
        permute=permute,
    )

    # Save results
    result_path = file_manager.results_dir(session_id, "combine") / "results.pkl"
    save_dataframe(df, result_path)

    log.info(f"Combine completed for session {session_id}: {len(df)} results")
    return result_path
