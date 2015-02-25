#! /usr/bin/env python
"""
Script for generation of input files for host to guest ITC titrations.
"""
from itctools.itctools import ureg, Quantity
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation

# Define solvents.
# TODO command line specification of density and name
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
buffer = Solvent('buffer', density=1.014 * ureg.gram / ureg.milliliter)

# The sample cell volume in microliters
cell_volume = 202.8

# Define compounds.

nguests = 8  # overnight from 5pm till 9am


host = Compound('host', molecular_weight=1162.9632 * ureg.gram / ureg.mole, purity=0.7133)
guest_molecular_weights = [
    209.12,
    123.62,
    153.65,
    189.13,
    187.11,
    151.63,
    135.64,
    149.66,
    163.69,
    238.59,
    147.65,
    189.73,
    173.68,
    203.71]
guest_compound_Ka = Quantity([21024287.6408,
                              13556262.0311,
                              81495.6444611,
                              1788684.70709,
                              2153596.60855,
                              769185.744612,
                              29967627.9188,
                              594389946.514,
                              2372114592.34,
                              683472.220385,
                              164811515.64,
                              6869559660.36,
                              28356538311.4,
                              396415131.021],
                             ureg.liter / ureg.mole)  # see dGtoKa.py
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


host_solution = SimpleSolution(
    compound=host,
    compound_mass=33.490 *
    ureg.milligram,
    solvent=buffer,
    solvent_mass=10.1151 *
    ureg.gram,
    location=PipettingLocation(
        source_plate.RackLabel,
        source_plate.RackType,
        1))
guest_solutions = list()

# Dispensed by quantos
guest_compound_masses = Quantity([2.210,
                                  1.400,
                                  1.705,
                                  1.945,
                                  1.975,
                                  1.635,
                                  1.700,
                                  1.640,
                                 ],
                                 ureg.milligram)
# Dispensed by quantos
guest_solvent_masses = Quantity([11.0478,
                                 11.6656,
                                 11.3653,
                                 10.8034,
                                 10.9705,
                                 10.8984,
                                 13.0758,
                                 10.9321,
                                ],
                                ureg.gram)

for guest_index in range(nguests):
    guest_solutions.append(
        SimpleSolution(
            compound=guests[guest_index],
            compound_mass=guest_compound_masses[guest_index],
            solvent=buffer,
            solvent_mass=guest_solvent_masses[guest_index],
            location=PipettingLocation(
                source_plate.RackLabel,
                source_plate.RackType,
                2 + guest_index)))

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
control_protocol.export_inj_file()
# Protocol for 1:1 binding analyis
blank_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Chodera Load Cell Without Cleaning Cell After.setup',
    itc_method='ChoderaHostGuest.inj',
    analysis_method='Onesite',
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
itc_experiment_set = ITCExperimentSet(name='SAMPL4-CB7 host-guest experiments')
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
optimal_rm = list()

for guest_index in range(nguests):

    # We need to store the experiments before adding them to the set
    host_guest_experiments = list()
    buff_guest_experiments = list()

    # Scaling factors per replicate
    factors = list()
    # Define host into guest experiments.
    for replicate in range(1):
        name = 'host into %s' % guests[guest_index].name
        experiment = ITCHeuristicExperiment(
            name=name,
            syringe_source=host_solution,
            cell_source=guest_solutions[guest_index],
            protocol=binding_protocol,
            cell_volume=cell_volume,
            cell_concentration=0.3 * (ureg.millimole / ureg.liter) * cell_scaling,
            buffer_source=buffer_trough)
        # optimize the syringe_concentration using heuristic equations and known binding constants
        optimal_rm.append(experiment.heuristic_syringe(
            guest_compound_Ka[guest_index],
            strict=False))
        # rescale if syringe > stock. Store factor.
        factors.append(experiment.rescale())
        host_guest_experiments.append(experiment)

    # Define buffer into guest experiments.
    for replicate in range(1):
        name = 'buffer into %s' % guests[guest_index].name
        experiment = ITCHeuristicExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=guest_solutions[guest_index],
            protocol=blank_protocol,
            cell_volume=cell_volume,
            cell_concentration=0.3 *
            ureg.millimole /ureg.liter,
            buffer_source=buffer_trough)
        # rescale to match host into guest experiment concentrations.
        experiment.rescale(tfactor=factors[replicate])
        buff_guest_experiments.append(experiment)

    # Add buffer to guest experiment(s) to set
    for buff_guest_experiment in buff_guest_experiments:
        itc_experiment_set.addExperiment(buff_guest_experiment)

    # Add host to guest experiment(s) to set
    for host_guest_experiment in host_guest_experiments:
        itc_experiment_set.addExperiment(host_guest_experiment)
        host_guest_experiment.simulate(guest_compound_Ka[guest_index], macromol_titrant=True,
                        filename='%s-simulation.png' % guests[guest_index].name)


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
excel_filename = 'host-guest-itc.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
