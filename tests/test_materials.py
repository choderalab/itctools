import unittest
from itctools import materials
from simtk import unit

class TestMaterials(unittest.TestCase):
    """Validation of the materials submodule"""

    def _check_solvent(self, args):
        """Validate the structure of a single or two argument Solvent."""
        defaults = ['', None]
        water = materials.Solvent(*args)
         # Pad the list with elements from defaults till it is the same length
        args.extend(defaults[len(args):])
        self.assertListEqual([water.name, water.density], args)

    def test_solvent(self):
        """Ensure Solvent has correct parameter assignment"""
        arguments = [['water'], ['water', 0.9970479 * unit.grams / unit.milliliter]]
        for args in arguments:
            self._check_solvent(args)

    def _check_compound(self, args):
        """Validate structure of single, two or three arguments"""
        defaults = ['', None, 1.00]
        compound = materials.Compound(*args)
        # Pad the list with elements fom defaults till it is the same length
        args.extend(defaults[len(args):])
        self.assertListEqual([compound.name,compound.molecular_weight,compound.purity, ], args)

    def test_compound(self):
        """Ensure Compound has correct parameter assignment"""
        arguments = [['nacl'],
                     ['imatinib mesylate', 589.7 * unit.grams/unit.mole],
                     ['compound1', 209.12 * unit.grams/unit.mole, 0.975]
        ]

        for args in arguments:
            self._check_compound(args)