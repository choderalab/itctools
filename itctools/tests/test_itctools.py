import unittest
from itctools import itctools as itc


class TestPermutation(unittest.TestCase):

    """Validation of the itctools.permutation_without_replacement function."""

    def _check_uniques(self, n, seq):
        """A set can't hold multiple copies of the same value.
        If the list only has uniques, it has the same length as the set.
        """
        compositions = itc.permutation_with_replacement(n, seq)
        self.assertEqual(len(compositions), len(set(compositions)))

    def test_permutation_with_replacement_uniques(self):
        """Ensure that permutation_with_replacement only generates unique sequences (when given unique inputs)."""
        seq = ['So test', 'Much robust', 'Very unit', 'Wow']
        for x in range(0, 7):
            self._check_uniques(x, seq)

    def _check_length(self, n, seq):
        """Makes sure that every sequence has the correct length.
        """
        compositions = itc.permutation_with_replacement(n, seq)
        for comp in compositions:
            self.assertEqual(len(comp), n)

    def test_permutation_with_replacement_length(self):
        """Ensure that permutation_with_replacement generates correct length."""
        seq = ['Did', 'you', 'mean', '"dog"?']
        for x in range(0, 7):
            self._check_length(x, seq)

    def test_ureg(self):
        """Ensure that the units are working"""
        quantity = itc.Quantity(140.0, 'micromolar')
        self.assertAlmostEqual(quantity / itc.ureg.standard_concentration, 1.4E-4, places=7)


class TestRm(unittest.TestCase):
    """Validation of the itctools.compute_rm function"""

    def test_Rm_over_1(self):
        """Ensure that Rm return values greater than 1 with strange pint units."""
        for value in [1000, 1000000., 5000000, 10000000]:
            c = itc.ureg.Quantity(value, 'millimole/mole')
            rm = itc.compute_rm(c)
            # if this fails, we are using weird "milli" units again
            self.assertGreaterEqual(itc.compute_rm(c), 1.0)

if __name__ == '__main__':
    unittest.main()
