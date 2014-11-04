import unittest
from itctools import materials
from simtk import unit


class TestMaterials(unittest.TestCase):

    """Validation of the materials submodule"""

    def _check_solvent(self, args):
        """Validate the structure of a single or two argument Solvent."""
        # [name, density=None]
        defaults = ['', None]
        water = materials.Solvent(*args)
        # Pad the list with defaults to account for missing arguments
        args.extend(defaults[len(args):])
        self.assertListEqual([water.name, water.density], args)

    def test_solvent(self):
        """Ensure Solvent has correct parameter assignment"""
        arguments = [
            ['water'],
            ['water', 0.9970479 * unit.grams / unit.milliliter],
            ]
        for args in arguments:
            self._check_solvent(args)

    def _check_compound(self, args):
        """Validate structure of single, two or three arguments"""
        # [name, molecular_weight=None, purity=1.0]
        defaults = ['', None, 1.00]
        compound = materials.Compound(*args)
        # Pad the list with defaults to account for missing arguments
        args.extend(defaults[len(args):])
        self.assertListEqual(
            [compound.name, compound.molecular_weight, compound.purity, ], args)

    def test_compound(self):
        """Ensure Compound has correct parameter assignment"""
        arguments = [['nacl'],
                     ['imatinib mesylate', 589.7 * unit.grams / unit.mole],
                     ['compound1', 209.12 * unit.grams / unit.mole, 0.975]
                     ]

        for args in arguments:
            self._check_compound(args)

    def _check_pureliquid(self, args):
        """Validate the structure of a triple or quadruple argument PureLiquid"""
        # [name, density, molecular_weight, purity=1.0]
        defaults = ['','','', 1.0]
        pureliquid = materials.PureLiquid(*args)
         # Pad the list with defaults to account for missing arguments
        args.extend(defaults[len(args):])
        self.assertListEqual(
            [pureliquid.name, pureliquid.density, pureliquid.molecular_weight, pureliquid.purity, ], args)

    def test_pure_liquid(self):
        """Ensure PureLiquid has correct parameter assignment"""
        arguments = [['water', 0.9970479 * unit.grams / unit.milliliter, 18.01528 * unit.gram / unit.mole],
                     ['ethanol', 0.789 * unit.grams / unit.milliliter, 46.07 * unit.gram / unit.mole, 99.8 / 100.]
                     ]
        for args in arguments:
            self._check_pureliquid(args)
