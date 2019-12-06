#! /usr/bin/env python
"""
Script for generation of input files for PHIPA ITC.

Ligand information is provided via a YAML file.
"""

from itctools.itctools import ureg, Quantity
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation
import pprint
import copy
import pandas as pd

# Read CSV file containing information about compounds
solubilities_df = pd.read_csv('swissadme_logS_53-compounds.csv')
solubilities_df.set_index('ID', inplace=True)
purities_df = pd.read_csv('compound_list.csv')
purities_df.set_index('Identifier', inplace=True)

# Compound IDs from CSV
compound_ids = [1, 2]

# The sample cell volume in microliters
cell_volume = 202.8

# Define solvents.
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
buffer = Solvent('buffer', density=1.014 * ureg.gram / ureg.milliliter) # TODO is our density the same as the HOST-GUEST buffer?

# Define compounds.
receptor = Compound('receptor', molecular_weight=29246.0 * (ureg.gram / ureg.mole), purity=1.)
ligands = list()
ligand_kas = list()
for compound_id in compound_ids:
    molecular_weight = solubilities_df.loc[compound_id]['MW'] * (ureg.gram / ureg.mole)
    purity = purities_df.loc[compound_id]['Purity'] / 100.0
    solubility = solubilities_df.loc[compound_id]['ESOL Solubility (mg/ml)'] * (ureg.milligram / ureg.milliliter)
    compound = Compound(f'{compound_id}', molecular_weight=molecular_weight, purity=purity, solubility=solubility)
    ligands.append(compound)
    Ka = 1.203e6 / ureg.molar
    ligand_kas.append(Ka)
    print(f"Compound {compound_id:8} : MW {(molecular_weight/(ureg.gram / ureg.mole)):8.3f} g/mol,  purity {purity:.3f}, solubility = {(solubility.to('milligrams/milliliter')):8.3f}, Kd = {((1/Ka).to('micromolar')):8.3f}")

# Define troughs on the instrument
water_trough = Labware(RackLabel='Water', RackType='Trough 100ml')
buffer_trough = Labware(RackLabel='Buffer', RackType='Trough 100ml')

# Define source labware.
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')
protein_source_plate = Labware(RackLabel='ProteinSourcePlate', RackType='ITC Plate')

# Define source solutions in the vial holder
protein_mass = (2.486 * ureg.milligram / ureg.milliliter) * (1.3 * ureg.milliliter)
solvent_mass = 1.3 * ureg.milliliters * buffer.density
receptor_solution = SimpleSolution(compound=receptor, compound_mass=protein_mass, solvent=buffer, solvent_mass=solvent_mass, location=PipettingLocation(
    protein_source_plate.RackLabel,
    protein_source_plate.RackType,
    1)) # Well A1 of ITC plate

ligand_solutions = list()
for ligand in ligands:
    ligand_solution = SimpleSolution(compound=ligand, compound_mass=10.0 * ureg.milligram, solvent=buffer, solvent_mass=10.0 * ureg.gram, location=PipettingLocation(
        source_plate.RackLabel,
        source_plate.RackType,
        1)) # Well A1 of vial holder
    ligand_solutions.append(ligand_solution)

# Define ITC protocol.

# CAII cell concentrations to evaluate
cell_concentrations = [0.010 * ureg.millimolar, 0.020 * ureg.millimolar, 0.040 * ureg.millimolar]

# Protocol for 'control' titrations (water-water, buffer-buffer,
# titrations into buffer, etc.)

control_protocol = ITCProtocol(
    'control_protocol',
    sample_prep_method='Plates Clean.setup', # thorough cleaning
    itc_method='ChoderaWaterWater.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )

# Protocol for 1:1 binding analysis
blank_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Chodera Load Cell Without Cleaning Cell After.setup',
    itc_method='ChoderaHSABlank.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )


binding_protocol = ITCProtocol(
    '1:1 binding protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='ChoderaHSA.inj',
    analysis_method='Onesite',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )
# Protocol for cleaning protocol
cleaning_protocol = ITCProtocol(
    'cleaning protocol',
    sample_prep_method='Plates Clean.setup', # thorough cleaning
    itc_method='water5inj.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=5 * [dict(volume_inj=7.5, duration_inj=15, spacing=150, filter_period=5)],
    )

# Define ITC Experiment.

# use specified protocol by default
itc_experiment_set = ITCExperimentSet(name='Receptor-ligand experiments')
# Add available plates for experiments.
destination_plate = Labware(
    RackLabel='DestinationPlate',
    RackType='ITC Plate')
itc_experiment_set.addDestinationPlate(destination_plate)

nreplicates = 2  # number of replicates of each experiment

# Add cleaning experiment.
name = 'initial cleaning water titration'
itc_experiment_set.addExperiment(
    ITCExperiment(
        name=name,
        syringe_source=water_trough,
        cell_source=water_trough,
        protocol=cleaning_protocol,
        cell_volume=cell_volume))

# Add water control titrations.
for replicate in range(1):
    name = 'water into water %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=buffer_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# buffer into receptor at all cell concentrations
for cell_concentration in cell_concentrations:
    for replicate in range(1):
        name = 'buffer into receptor %d' % (replicate + 1)
        itc_experiment_set.addExperiment(
            ITCExperiment(
                name=name,
                syringe_source=buffer_trough,
                cell_source=receptor_solution,
                protocol=control_protocol,
                cell_concentration=cell_concentration,
                buffer_source=buffer_trough,
                cell_volume=cell_volume))

# ligands/receptor
# scale cell concentration to fix necessary syringe concentrations

for ligand, ligand_solution, ligand_ka in zip(ligands, ligand_solutions, ligand_kas):

    # We need to store the experiments before adding them to the set
    ligand_protein_experiments = list()
    ligand_buffer_experiments = list()

    # Scaling factors per replicate
    factors = list()

    # Define ligand to protein experiments
    for cell_concentration in cell_concentrations:
        for replicate in range(nreplicates):
            name = '%s into receptor %d' % (ligand.name, replicate + 1 )
            experiment = ITCHeuristicExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=receptor_solution,
                protocol=binding_protocol,
                cell_concentration=cell_concentration,
                buffer_source=buffer_trough,
                cell_volume=cell_volume)
            # optimize the syringe_concentration using heuristic equations and known binding constants
            # TODO extract m, v and V0 from protocol somehow?

            # Warning, you're possibly not getting the setup you want. Consider not using the Heuristic Experiment
            experiment.heuristic_syringe(ligand_ka, 10, strict=False)
            # rescale if syringe > stock. Store factor.
            factors.append(experiment.rescale())
            ligand_protein_experiments.append(experiment)

    # Define corresponding ligand into buffer experiments
    for experiment in ligand_protein_experiments:
        experiment = copy.deepcopy(experiment)
        experiment.cell_source = buffer_trough
        experiment.cell_concentration *= 0
        ligand_buffer_experiments.append(experiment)

    # TODO, since we are changing ligands, we'd have to wash the syringe.
    # Add ligand to protein experiment(s) to set
    for ligand_protein_experiment in ligand_protein_experiments:
        itc_experiment_set.addExperiment(ligand_protein_experiment)
        # pprint.pprint(ligand_protein_experiment.__dict__)

    # Add ligand_to_buffer experiment(s) to set
    for ligand_buffer_experiment in ligand_buffer_experiments:
        itc_experiment_set.addExperiment(ligand_buffer_experiment)
        # pprint.pprint(ligand_buffer_experiment.__dict__)


# Add cleaning experiment.
name = 'final cleaning water titration'
itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=cleaning_protocol, cell_volume=cell_volume) )

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
            cell_volume=cell_volume))

# Check that the experiment can be carried out using available solutions and plates.
# Also generate Tecan EVO worklists
import sys
itc_experiment_set.validate(print_volumes=True, omit_zeroes=True, human_readable_log=sys.stdout)

# For convenience, concentrations
for ligand_solution in ligand_solutions:
    print("%12s %.4f mM" % (ligand_solution.name, ligand_solution.concentration / ureg.millimolar ))
    print("%12s %.4f mM" % ('receptor', receptor_solution.concentration  / ureg.millimolar ))


# Write Tecan EVO pipetting operations for both liquid and air LiHas.
worklist_prefix = 'experiments'
itc_experiment_set.writeTecanWorklist(worklist_prefix)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'experiments.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
