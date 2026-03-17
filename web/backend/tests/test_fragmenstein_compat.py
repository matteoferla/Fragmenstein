"""Compatibility tests — verify our backend assumptions against the actual Fragmenstein library.

These tests catch breaking changes in the upstream Fragmenstein API that would
silently break the web UI. They test the library contracts we depend on, not our code.
"""

import pytest
from rdkit import Chem


# ── Public API imports ────────────────────────────────────────────────

def test_public_imports():
    """Verify all classes we import from fragmenstein still exist."""
    from fragmenstein import Victor, Laboratory, Wictor, Monster
    assert Victor is not None
    assert Laboratory is not None
    assert Wictor is not None
    assert Monster is not None


def test_optional_imports():
    """Optional classes should either import or raise ImportError cleanly."""
    try:
        from fragmenstein import Quicktor
        assert Quicktor is not None
    except ImportError:
        pass  # acceptable

    try:
        from fragmenstein import OpenVictor
        assert OpenVictor is not None
    except ImportError:
        pass  # acceptable


# ── Victor API contract ──────────────────────────────────────────────

def test_victor_class_attributes():
    """Verify Victor has the class attributes we depend on."""
    from fragmenstein import Victor
    # Pre-existing attributes
    assert hasattr(Victor, "work_path")
    assert hasattr(Victor, "monster_throw_on_discard")
    assert hasattr(Victor, "quick_reanimation")
    assert hasattr(Victor, "error_to_catch")
    assert hasattr(Victor, "uses_pyrosetta")
    # monster_joining_cutoff is set dynamically — verify it's settable
    Victor.monster_joining_cutoff = 5.0
    assert Victor.monster_joining_cutoff == 5.0


def test_victor_extract_mol_signature():
    """Verify Victor.extract_mol accepts the parameters we pass."""
    from fragmenstein import Victor
    import inspect
    sig = inspect.signature(Victor.extract_mol)
    params = list(sig.parameters.keys())
    assert "name" in params
    assert "filepath" in params or "block" in params
    assert "smiles" in params
    assert "ligand_resn" in params
    assert "proximityBonding" in params


def test_victor_combine_signature():
    """Verify Victor.combine accepts warhead_harmonisation."""
    from fragmenstein import Victor
    import inspect
    sig = inspect.signature(Victor.combine)
    params = list(sig.parameters.keys())
    assert "warhead_harmonisation" in params


def test_victor_place_signature():
    """Verify Victor.place accepts merging_mode."""
    from fragmenstein import Victor
    import inspect
    sig = inspect.signature(Victor.place)
    params = list(sig.parameters.keys())
    assert "merging_mode" in params


def test_victor_summarize_keys():
    """Verify Victor.summarize returns the expected keys by checking a real run."""
    from fragmenstein import Wictor
    from fragmenstein.demo import MPro

    template = MPro.get_template()
    hits = [MPro.get_mol("x0387"), MPro.get_mol("x0434")]

    Wictor.work_path = "/tmp/frag_compat_test"
    Wictor.error_to_catch = Exception

    v = Wictor(hits=hits, pdb_block=template, ligand_resn="LIG", covalent_resi="145A")
    v.combine()
    summary = v.summarize()

    expected_keys = {"name", "smiles", "error", "mode", "∆∆G", "∆G_bound", "∆G_unbound",
                     "comRMSD", "N_constrained_atoms", "N_unconstrained_atoms", "runtime",
                     "regarded", "disregarded"}
    missing = expected_keys - set(summary.keys())
    assert not missing, f"Victor.summarize() missing keys: {missing}"


# ── Laboratory API contract ──────────────────────────────────────────

def test_laboratory_constructor_signature():
    """Verify Laboratory.__init__ accepts the params we pass."""
    from fragmenstein import Laboratory
    import inspect
    sig = inspect.signature(Laboratory.__init__)
    params = list(sig.parameters.keys())
    assert "pdbblock" in params
    assert "covalent_resi" in params
    assert "run_plip" in params


def test_laboratory_combine_signature():
    """Verify Laboratory.combine accepts what we pass."""
    from fragmenstein import Laboratory
    import inspect
    sig = inspect.signature(Laboratory.combine)
    params = list(sig.parameters.keys())
    assert "mols" in params
    assert "permute" in params
    assert "combination_size" in params


def test_laboratory_place_signature():
    """Verify Laboratory.place accepts queries."""
    from fragmenstein import Laboratory
    import inspect
    sig = inspect.signature(Laboratory.place)
    params = list(sig.parameters.keys())
    assert "queries" in params


def test_laboratory_category_labels():
    """Verify outcome categories haven't changed."""
    from fragmenstein import Laboratory
    assert hasattr(Laboratory, "category_labels")
    labels = Laboratory.category_labels
    # These are the categories we use in the frontend
    for expected in ["acceptable", "crashed", "timeout"]:
        assert expected in labels, f"Missing category: {expected}"


def test_laboratory_victor_swappable():
    """Verify we can swap Laboratory.Victor class attribute."""
    from fragmenstein import Laboratory, Wictor
    original = Laboratory.Victor
    Laboratory.Victor = Wictor
    assert Laboratory.Victor is Wictor
    Laboratory.Victor = original


# ── Laboratory subprocess methods ────────────────────────────────────

def test_combine_subprocess_exists():
    """Verify combine_subprocess is an instance method we can override."""
    from fragmenstein import Laboratory
    assert hasattr(Laboratory, "combine_subprocess")
    assert callable(Laboratory.combine_subprocess)


def test_place_subprocess_exists():
    """Verify place_subprocess is an instance method we can override."""
    from fragmenstein import Laboratory
    assert hasattr(Laboratory, "place_subprocess")
    assert callable(Laboratory.place_subprocess)


# ── WebLaboratory subclass ───────────────────────────────────────────

def test_weblaboratory_subclass():
    """Verify our WebLaboratory subclass works with the current Library."""
    from app.services.combine_service import WebLaboratory
    from fragmenstein import Laboratory

    assert issubclass(WebLaboratory, Laboratory)
    assert hasattr(WebLaboratory, "_warhead_harmonisation")
    assert hasattr(WebLaboratory, "_merging_mode")


# ── Demo data ────────────────────────────────────────────────────────

def test_demo_mpro_data():
    """Verify MPro demo data is accessible."""
    from fragmenstein.demo import MPro

    template = MPro.get_template()
    assert isinstance(template, str)
    assert "ATOM" in template

    hits = MPro.get_n_filtered_mols(3)
    assert len(hits) == 3
    for mol in hits:
        assert isinstance(mol, Chem.Mol)
        assert mol.GetNumAtoms() > 0


def test_demo_mac1_data():
    """Verify Mac1 demo data is accessible."""
    from fragmenstein.demo import Mac1

    template = Mac1.get_template()
    assert isinstance(template, str)
    assert "ATOM" in template


# ── Monster standalone ───────────────────────────────────────────────

def test_monster_combine_api():
    """Verify Monster.combine works standalone (no protein)."""
    from fragmenstein import Monster
    from fragmenstein.demo import MPro

    hits = [MPro.get_mol("x0387"), MPro.get_mol("x0434")]
    monster = Monster(hits=hits)
    monster.combine()

    assert monster.positioned_mol is not None
    assert isinstance(monster.positioned_mol, Chem.Mol)
    assert monster.positioned_mol.GetNumAtoms() > 0


# ── Result DataFrame schema ──────────────────────────────────────────

def test_wictor_combine_result_columns():
    """Verify a Wictor combine produces the DataFrame columns we serialize."""
    from fragmenstein import Laboratory, Wictor
    from fragmenstein.demo import MPro

    template = MPro.get_template()
    hits = [MPro.get_mol("x0387"), MPro.get_mol("x0434")]

    Laboratory.Victor = Wictor
    Wictor.work_path = "/tmp/frag_compat_test"
    Wictor.error_to_catch = Exception

    lab = Laboratory(pdbblock=template, covalent_resi="145A")
    df = lab.combine(hits, n_cores=1, timeout=60)

    # Columns our result_serializer depends on
    expected = {"name", "smiles", "error", "mode", "∆∆G", "outcome"}
    present = set(df.columns)
    missing = expected - present
    assert not missing, f"Combine result DataFrame missing columns: {missing}"

    # simple_smiles is generated by combine
    assert "simple_smiles" in df.columns, "simple_smiles column missing — SmallWorld search would fail"
