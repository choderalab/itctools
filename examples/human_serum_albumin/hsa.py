#! /usr/bin/env python
"""
Script for generation of input files for Aspirin and Naproxen binding to HSA by ITC.
"""

from simtk.unit import *
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation

# Define solvents.
water = Solvent('water', density=0.9970479 * grams / milliliter)
buffer = Solvent('buffer', density=1.014 * grams / milliliter) # TODO is our density the same as the HOST-GUEST buffer?

# Define compounds.
hsa = Compound('HumanSerumAlbumin', molecular_weight=65000 * daltons, purity=.95)
aspirin = Compound('AcetylsalicylicAcid', molecular_weight=180.15742 * daltons, purity=.99)
naproxen = Compound('Naproxen', molecular_weight=230.3 * daltons, purity=.98)


# Define troughs on the instrument.

water_trough = Labware(RackLabel='Water', RackType='Trough 100ml')
buffer_trough = Labware(RackLabel='Buffer', RackType='Trough 100ml')

# Define source labware.
#source_plate = Labware(RackLabel='SourcePlate', RackType='12WellVialHolder')
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')

# Define source solutions on the deck.one
# TODO : Define solutions once prepared with the Quantos



# Define ITC protocol.

# Protocol for 'control' titrations (water-water, buffer-buffer,
# titrations into buffer, etc.)
control_protocol = ITCProtocol(
    'control protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaWaterWater.inj',
    analysis_method='Control')
# Protocol for 1:1 binding analyis
blank_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Chodera Load Cell Without Cleaning Cell After.setup',
    itc_method='ChoderaHostGuest.inj',  #  TODO Define new protocol?
    analysis_method='Onesite')
binding_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaHostGuest.inj',  #  TODO Define new protocol?
    analysis_method='Onesite')
# Protocol for cleaning protocol
cleaning_protocol = ITCProtocol(
    'cleaning protocol',
    sample_prep_method='Plates Clean.setup',
    itc_method='water5inj.inj',
    analysis_method='Control')

# Define ITC Experiment.

# use specified protocol by default
itc_experiment_set = ITCExperimentSet(name='Human Serum Albumin experiments')
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
        protocol=cleaning_protocol))

# Add water control titrations.
for replicate in range(1):
    name = 'water into water %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol))

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=buffer_trough,
            protocol=control_protocol))

# drug into buffer.
for replicate in range(1):
    name = 'host into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=XXXXX, #TODO define drug solutions
            cell_source=buffer_trough,
            protocol=binding_protocol))

# drugs/HSA
# scale cell concentration to fix necessary syringe concentrations

#TODO, make this loop run over the drugs
cell_scaling = 1.
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
            cell_concentration=0.2 *
            millimolar *
            cell_scaling,
            buffer_source=buffer_trough)
        # optimize the syringe_concentration using heuristic equations and known binding constants
        # TODO extract m, v and V0 from protocol somehow?
        experiment.heuristic_syringe(
            guest_compound_Ka[guest_index],
            10,
            3. *
            microliters,
            202.8 *
            microliters)
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
            cell_concentration=0.2 *
            millimolar,
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
            protocol=control_protocol))

# Check that the experiment can be carried out using available solutions
# and plates.

itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)

# For convenience, concentrations
for g, guest in enumerate(guest_solutions, start=1):
    print("guest%02d" % g, guest.concentration.in_units_of(millimolar))

print("host", host_solution.concentration.in_units_of(millimolar))


# Write Tecan EVO pipetting operations.
worklist_filename = 'hsa.gwl'
itc_experiment_set.writeTecanWorklist(worklist_filename)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'hsa.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
