#! /usr/bin/env python
"""
Script for High Dynamic Range itc experiment.

"""

import numpy as np
from itctools.itctools import ureg, Quantity
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation



#######################################
# Below are the params needed to specify.
#######################################


# Define solvents.
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
buffer = Solvent('buffer', density=1.014 * ureg.gram / ureg.milliliter)

# The sample cell volume in microliters
cell_volume = 202.8

# Define compounds.

# each experiment takes about 50 mins, cleaning takes 2.5
nguests = 1
assert(nguests==1)

host = Compound('host', molecular_weight=1134.987 * ureg.gram / ureg.mole, purity=0.98)
guest_molecular_weights = [123.11]


# If HDR is not activated, there is only one value to guest_compound_Ka.
# guest_compound_Ka = Quantity([np.exp(32)], ureg.liter / ureg.mole)  # see dGtoKa.py



# c = [M]_0 * Ka
# R_m = 6.4/c^0.2 + 13/c
# Hence the scalling of c is eqivalent to the scalling of the association constant.

guest_compound_Ka_guess = Quantity([np.exp(6), np.exp(8), np.exp(10), np.exp(12), np.exp(14)], ureg.liter / ureg.mole)

# We only allow one guest to be added into the experiment, for now.
guests = [
    Compound(
        name='guest%02d' %
        (guest_index +
         1),
        molecular_weight=guest_molecular_weights[guest_index] *
        ureg.gram / ureg.mole,
        purity=0.975) for guest_index in range(nguests)]



# Define troughs on the instrument.

water_trough = Labware(RackLabel='Water', RackType='Trough 100ml')
buffer_trough = Labware(RackLabel='Buffer', RackType='Trough 100ml')

# Define source labware.
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')

# Define source solutions on the deck.one
# TODO : Use actual compound and solvent masses.
# NOTE: Host solution is diluted by 10x.
# NOTE: it might not be a good idea to further dilute.
#       since the host solution is already very dilute because of the low solubility.


host_solution = SimpleSolution(
    compound=host,
    compound_mass=11.54 *
    ureg.milligram,
    solvent=buffer,
    solvent_mass=10.2002 *
    ureg.gram,
    location=PipettingLocation(
        source_plate.RackLabel,
        source_plate.RackType,
        1))
guest_solutions = list()

# Dispensed by quantos
guest_compound_masses = Quantity([1.38],
                                 ureg.milligram)
# Dispensed by quantos
guest_solvent_masses = Quantity([11.4985],
                                ureg.gram)


for ka_idx, ka in enumerate(list(guest_compound_Ka_guess)):
    guest_solutions.append(
        SimpleSolution(
            compound=guests[0],
            compound_mass=guest_compound_masses[0],
            solvent=buffer,
            solvent_mass=guest_solvent_masses[0],
            location=PipettingLocation(
                source_plate.RackLabel,
                source_plate.RackType,
                2 + ka_idx)))


# Have one solution just to react with buffer.



#######################################
# Above are the params needed to specify.
#######################################

# Define ITC protocol.

# Protocol for 'control' titrations (water-water, buffer-buffer,
# titrations into buffer, etc.)
control_protocol = ITCProtocol(
    'control_protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaWaterWater.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )

# Protocol for 1:1 binding analyis
blank_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaHostGuest.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )
binding_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaHostGuest.inj',
    analysis_method='Onesite',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )
# Protocol for cleaning protocol
cleaning_protocol = ITCProtocol(
    'cleaning protocol',
    sample_prep_method='Plates Clean.setup',
    itc_method='water5inj.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=5 * [dict(volume_inj=7.5, duration_inj=15, spacing=150, filter_period=5)],
    )


# Define ITC Experiment.

# use specified protocol by default
itc_experiment_set = ITCExperimentSet(name='hdr_experiments', hdr = True)
# Add available plates for experiments.
itc_experiment_set.addDestinationPlate(
    Labware(
        RackLabel='DestinationPlate',
        RackType='ITC Plate'))
itc_experiment_set.addDestinationPlate(
    Labware(
        RackLabel='DestinationPlate2',
        RackType='ITC Plate'))

nreplicates = 1  # number of replicates of each experiment

# Add cleaning experiment.
name = 'initial cleaning water titration'
itc_experiment_set.addExperiment(
    ITCExperiment(
        name=name,
        syringe_source=water_trough,
        cell_source=water_trough,
        protocol=cleaning_protocol,
        cell_volume=cell_volume,))

# Add water control titrations.
for replicate in range(1):
    name = 'water into water %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume,))

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=buffer_trough,
            protocol=control_protocol,
            cell_volume=cell_volume,))

# Host into buffer.
for replicate in range(1):
    name = 'host into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=host_solution,
            cell_source=buffer_trough,
            protocol=binding_protocol,
            cell_volume=cell_volume,))

# Host/guests.
# scale cell concentration to fix necessary syringe concentrations
cell_scaling = 1.
cell_concentration = 0.3 * ureg.millimole / ureg.liter * cell_scaling
optimal_rm = list()

for kda_idx, kda in enumerate(list(guest_compound_Ka_guess)):

    # We need to store the experiments before adding them to the set
    host_guest_experiments = list()
    buff_guest_experiments = list()

    # Scaling factors per replicate
    factors = list()
    # Define host into guest experiments.
    for replicate in range(1):
        name = 'host into guest %s' % kda_idx
        experiment = ITCHeuristicExperiment(
            name=name,
            syringe_source=host_solution,
            cell_source=guest_solutions[kda_idx],
            protocol=binding_protocol,
            cell_volume=cell_volume,
            cell_concentration=cell_concentration,
            buffer_source=buffer_trough)
        # optimize the syringe_concentration using heuristic equations and known binding constants
        optimal_rm.append(experiment.heuristic_syringe(
            guest_compound_Ka_guess[kda_idx],
            strict=False))
        # rescale if syringe > stock. Store factor.
        factors.append(experiment.rescale())
        host_guest_experiments.append(experiment)
    # Add host to guest experiment(s) to set
    for host_guest_experiment in host_guest_experiments:
        itc_experiment_set.addExperiment(host_guest_experiment)
        host_guest_experiment.simulate(guest_compound_Ka_guess[kda_idx], macromol_titrant=True,
                        filename='simulation.png')



# Define buffer into guest experiments.
for replicate in range(1):
    name = 'buffer into guest'
    experiment = ITCHeuristicExperiment(
        name=name,
        syringe_source=buffer_trough,
        cell_source=guest_solutions[-1],
        protocol=blank_protocol,
        cell_volume=cell_volume,
        cell_concentration=cell_concentration,
        buffer_source=buffer_trough)
    # rescale to match host into guest experiment concentrations.
    experiment.rescale(tfactor=factors[replicate])
    buff_guest_experiments.append(experiment)

# Add buffer to guest experiment(s) to set
for buff_guest_experiment in buff_guest_experiments:
    itc_experiment_set.addExperiment(buff_guest_experiment)


# Add cleaning experiment.
#name = 'final cleaning water titration'
#itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=cleaning_protocol) )

# Water control titrations.
nfinal = 2
for replicate in range(nfinal):
    name = 'final water into water test %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume,))

# Check that the experiment can be carried out using available solutions
# and plates.

itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)

# Get effective Rm
actual_Rm = list()
for experiment in itc_experiment_set.experiments:
    try:
        guest_mol = ((experiment.cell_concentration * 202.8 * Quantity('microliter')).to('millimole'))
        host_mol = ((experiment.syringe_concentration * 30 * Quantity('microliter')).to('millimole'))
        Rm = host_mol/guest_mol
        actual_Rm.append(Rm)

    except (AttributeError, TypeError):
        continue

with open("host-guest-itc-Rm.txt", 'w') as ratio_file:
    ratio_file.write("#, Optimal, Actual\n")
    for idx,(opt,act) in enumerate(zip(optimal_rm, actual_Rm), start=1):
        ratio_file.write("%d, %.5f, %.5f\n" % (idx, opt, act))


# Write Tecan EVO pipetting operations.
worklist_filename = 'host-guest-itc.gwl'
itc_experiment_set.writeTecanWorklist(worklist_filename)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'host-guest-itc.csv'
itc_experiment_set.writeAutoITCCSV(excel_filename)
