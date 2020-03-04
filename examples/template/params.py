#! /usr/bin/env python

from template_PHIPA import generate_files
from itctools.itctools import ureg, Quantity

params = {
    "protein_volume": 0.5 * ureg.milliliter,
    "protein_concentration": 1.872 * ureg.milligram / ureg.milliliter,
    "buffer_masses": {'Compound 10' : 10.065 * ureg.gram},
    "compound_masses": {'Compound 10' : 15.05 * ureg.milligram},
    "cell_concentrations": [0.025 * ureg.millimolar],
    "number_of_itc_plates": 2
}
generate_files(**params)

