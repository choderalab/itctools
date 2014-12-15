import itertools
from pint import UnitRegistry

ureg = UnitRegistry()
Quantity = ureg.Quantity
ureg.define('molar = 1 * mole / liter = M')
ureg.define('standard_concentration = 1 M')

def permutation_with_replacement(n, seq):
    """
    Returns a list of all possible combinations of elements in seq, with length n.
    (Like permutation with replacement)
    """
    options = list()
    for p in itertools.product(seq, repeat=n):
        options.append(p)
    return options

