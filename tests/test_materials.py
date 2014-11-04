import unittest
from itctools import materials
from simtk import unit
from mock import Mock


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
            [compound.name, compound.molecular_weight, compound.purity, ],
            args)

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
        defaults = ['', '', '', 1.0]
        pureliquid = materials.PureLiquid(*args)
        # Pad the list with defaults to account for missing arguments
        args.extend(defaults[len(args):])
        self.assertListEqual([pureliquid.name,
                              pureliquid.density,
                              pureliquid.molecular_weight,
                              pureliquid.purity,
                              ],
                             args)

    def test_pure_liquid(self):
        """Ensure PureLiquid has correct parameter assignment"""
        arguments = [['water', 0.9970479 *
                      unit.grams /
                      unit.milliliter, 18.01528 *
                      unit.gram /
                      unit.mole], ['ethanol', 0.789 *
                                   unit.grams /
                                   unit.milliliter, 46.07 *
                                   unit.gram /
                                   unit.mole, 99.8 /
                                   100.]]
        for args in arguments:
            self._check_pureliquid(args)

    def test_simple_solution(self):
        """Ensure that SimpleSolution has correct parameter assignment"""

        # Using Mock classes to ensure that SimpleSolution works, independently
        # of Compound and Solvent
        mock_compound = Mock(
            name='imatinib',
            molecular_weight=589.7 *
            unit.grams /
            unit.mole,
            purity=1.0)
        compound_mass = 1.0 * unit.milligram
        mock_solvent = Mock(
            name='water', density=0.9970479 * unit.grams / unit.milliliter)
        solvent_mass = 10.0 * unit.grams
        mock_location = Mock(
            RackLabel='DestinationPlate', RackType='ITCPlate', Position=1)

        simple_solution = materials.SimpleSolution(
            mock_compound,
            compound_mass,
            mock_solvent,
            solvent_mass,
            mock_location)

        expected_structure = {
            'solvent': mock_solvent,
            'name': mock_compound.name,
            'compound_mass': compound_mass,
            'density': mock_solvent.density,  # Ideal solution approximation
            'volume': solvent_mass / mock_solvent.density,
            'location': mock_location,
            'compound': mock_compound,
            'solution_mass': compound_mass + solvent_mass,
            'compound_moles': compound_mass / mock_compound.molecular_weight * mock_compound.purity,
            'solvent_mass': solvent_mass,
            'concentration': (compound_mass / mock_compound.molecular_weight * mock_compound.purity) / (solvent_mass / mock_solvent.density),
        }
        try:  # python 3.2 +
            self.assertCountEqual(
                expected_structure, simple_solution.__dict__)
        except AttributeError:
            self.assertItemsEqual(
                expected_structure, simple_solution.__dict__)

    def test_simple_solution_solvent_compatibility(self):
        """Ensure that SimpleSolution is compatible with Solvent"""

        # Using Mock classes to ensure that SimpleSolution works, independently
        # of Compound
        mock_compound = Mock(
            name='imatinib',
            molecular_weight=589.7 *
            unit.grams /
            unit.mole,
            purity=1.0)
        compound_mass = 1.0 * unit.milligram
        solvent = materials.Solvent(
            'water', density=0.9970479 * unit.grams / unit.milliliter)
        solvent_mass = 10.0 * unit.grams
        mock_location = Mock(
            RackLabel='DestinationPlate', RackType='ITCPlate', Position=1)

        simple_solution = materials.SimpleSolution(
            mock_compound, compound_mass, solvent, solvent_mass, mock_location)

        expected_structure = {
            'solvent': solvent,
            'name': mock_compound.name,
            'compound_mass': compound_mass,
            'density': solvent.density,  # Ideal solution approximation
            'volume': solvent_mass / solvent.density,
            'location': mock_location,
            'compound': mock_compound,
            'solution_mass': compound_mass + solvent_mass,
            'compound_moles': compound_mass / mock_compound.molecular_weight * mock_compound.purity,
            'solvent_mass': solvent_mass,
            'concentration': (compound_mass / mock_compound.molecular_weight * mock_compound.purity) / (solvent_mass / solvent.density),
        }

        try:  # python 3.2 +
            self.assertCountEqual(
                expected_structure, simple_solution.__dict__)
        except AttributeError:
            self.assertItemsEqual(
                expected_structure, simple_solution.__dict__)

    def test_simple_solution_compound_compatibility(self):
        """Ensure that SimpleSolution is compatible with Compound"""

        # Using Mock classes to ensure that SimpleSolution works, independently
        # of Compound and Solvent
        compound = materials.Compound(
            'imatinib',
            molecular_weight=589.7 *
            unit.grams /
            unit.mole,
            purity=1.0)
        compound_mass = 1.0 * unit.milligram
        mock_solvent = Mock(
            name='water', density=0.9970479 * unit.grams / unit.milliliter)
        solvent_mass = 10.0 * unit.grams
        mock_location = Mock(
            RackLabel='DestinationPlate', RackType='ITCPlate', Position=1)

        simple_solution = materials.SimpleSolution(
            compound, compound_mass, mock_solvent, solvent_mass, mock_location)

        expected_structure = {
            'solvent': mock_solvent,
            'name': compound.name,
            'compound_mass': compound_mass,
            'density': mock_solvent.density,  # Ideal solution approximation
            'volume': solvent_mass / mock_solvent.density,
            'location': mock_location,
            'compound': compound,
            'solution_mass': compound_mass + solvent_mass,
            'compound_moles': compound_mass / compound.molecular_weight * compound.purity,
            'solvent_mass': solvent_mass,
            'concentration': (compound_mass / compound.molecular_weight * compound.purity) / (solvent_mass / mock_solvent.density),
        }

        try:  # python 3.2+
            self.assertCountEqual(
                expected_structure, simple_solution.__dict__)
        except AttributeError:
            self.assertItemsEqual(
                expected_structure, simple_solution.__dict__)
