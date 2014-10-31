import unittest
from itctools import materials
from simtk import unit

class TestMaterials(unittest.TestCase):
    """Validation of the materials submodule"""

    def _check_solvent(self, args):
        water = materials.Solvent(*args)
        if len(args) < 2:
            args.append(None)
        self.assertListEqual([water.name, water.density], args)

    def test_solvent(self):
        """Ensure Solvent has correct parameter assignment"""
        arguments = [['water'], ['water', 0.9970479 * unit.grams / unit.milliliter]]
        for args in arguments:
            self._check_solvent(args)
