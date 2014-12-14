import unittest
from itctools import itctools as itc


class TestITC(unittest.TestCase):

    """Validation of the itctools submodule."""

    def _check_uniques(self, n, seq):
        """A set can't hold multiple copies of the same value.
        If the list only has uniques, it has the same length as the set.
        """
        compositions = itc.permutation_with_replacement(n, seq)
        self.assertEqual(len(compositions), len(set(compositions)))

    def test_permutation_with_replacement_uniques(self):
        """Make sure that permutation_with_replacement only generates unique sequences (when given unique inputs)."""
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
        """Make sure that permutation_with_replacement generates correct length."""
        seq = ['Did', 'you', 'mean', '"dog"?']
        for x in range(0, 7):
            self._check_length(x, seq)

if __name__ == '__main__':
    unittest.main()
