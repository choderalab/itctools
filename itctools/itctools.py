import itertools
from pint import UnitRegistry
import numpy

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

@ureg.wraps(ureg.dimensionless, ureg.dimensionless, strict=False)
def compute_rm(c):
    rm = 6.4 / numpy.power(c, 0.2) + 13 / c
    return rm
