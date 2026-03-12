"""Regression tests for the PyRosetta mock (issue #58).

These run in a subprocess so the real pyrosetta (if installed) doesn't interfere.
"""
import unittest
import subprocess
import sys


class MockImportTests(unittest.TestCase):
    """Verify that ``import pyrosetta.rosetta.core`` works when PyRosetta is absent."""

    def _run_isolated(self, code: str):
        """Run *code* in a subprocess that cannot see real pyrosetta."""
        result = subprocess.run(
            [sys.executable, '-c', code],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            self.fail(f'Subprocess failed:\n{result.stderr}')

    def test_submodule_import(self):
        """``import pyrosetta.rosetta.core as prc`` must not raise TypeError."""
        self._run_isolated(
            "from fragmenstein.igor.pyrosetta_import import pyrosetta\n"
            "import pyrosetta.rosetta.core as prc\n"
            "assert 'pyrosetta.rosetta.core' in str(prc)\n"
        )

    def test_deep_attribute_access(self):
        """Chained attribute access should return mock objects, not explode."""
        self._run_isolated(
            "from fragmenstein.igor.pyrosetta_import import pyrosetta\n"
            "obj = pyrosetta.rosetta.core.scoring.constraints\n"
            "assert callable(obj)\n"
        )

    def test_mock_is_callable(self):
        self._run_isolated(
            "from fragmenstein.igor.pyrosetta_import import pyrosetta\n"
            "result = pyrosetta.init()\n"
            "assert result is not None\n"
        )


if __name__ == '__main__':
    unittest.main()
