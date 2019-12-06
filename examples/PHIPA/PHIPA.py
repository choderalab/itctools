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
compound_ids = [
    46, # acetyllysine
    #24, # displaces water network
    ]

nreplicates = 2  # number of replicates of each protein-ligand experiment

# The sample cell volume in microliters
cell_volume = 202.8

# Define solvents.
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
buffer = Solvent('buffer', density=1.014 * ureg.gram / ureg.milliliter) # TODO is our density the same as the HOST-GUEST buffer?

# Define receptor
receptor = Compound('receptor', molecular_weight=14117.0 * (ureg.gram / ureg.mole), purity=1.)

# Define compounds.
ligands = list()
ligand_kas = list()
print('COMPOUNDS')
for compound_id in compound_ids:
    molecular_weight = solubilities_df.loc[compound_id]['MW'] * (ureg.gram / ureg.mole)
    purity = purities_df.loc[compound_id]['Purity'] / 100.0
    solubility = solubilities_df.loc[compound_id]['ESOL Solubility (mg/ml)'] * (ureg.milligram / ureg.milliliter)
    compound = Compound(f'Compound {compound_id}', molecular_weight=molecular_weight, purity=purity, solubility=solubility)
    ligands.append(compound)
    Ka = 1.203e6 / ureg.molar
    ligand_kas.append(Ka)
    print(f"{compound.name:20} : MW {(molecular_weight/(ureg.gram / ureg.mole)):8.3f} g/mol,  purity {purity:.3f}, solubility = {(solubility.to('milligrams/milliliter')):8.3f}, Kd = {((1/Ka).to('micromolar')):8.3f}")
print('')

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

solubility_safety = 0.5 # fraction of solubility to prepare solution for
max_compound_mass = 10.0 * ureg.milligram

ligand_solutions = list()
for ligand in ligands:
    target_concentration = solubility_safety * ligand.solubility
    target_buffer_mass = 10.0 * ureg.gram
    target_buffer_volume = target_buffer_mass / buffer.density
    target_compound_mass = target_concentration * target_buffer_volume
    target_compound_mass = min(target_compound_mass, max_compound_mass)
    ligand_solution = SimpleSolution(compound=ligand, compound_mass=target_compound_mass, solvent=buffer, solvent_mass=target_buffer_mass, location=PipettingLocation(
        source_plate.RackLabel,
        source_plate.RackType,
        1)) # Well A1 of vial holder
    ligand_solutions.append(ligand_solution)

# For convenience, report concentrations
print('STOCK SOLUTION CONCENTRATIONS:')
print(f"{receptor_solution.name:20} : {receptor_solution.compound_mass.to('milligrams')} receptor in {receptor_solution.solvent_mass.to('grams')} solvent : {receptor_solution.concentration.to('millimolar'):30}")
for ligand_solution in ligand_solutions:
    print(f"{ligand_solution.name:20} : {ligand_solution.compound_mass.to('milligrams')} compound in {ligand_solution.solvent_mass.to('grams')} solvent : {ligand_solution.concentration.to('millimolar'):30}")
print('')

# Define ITC protocol.

# Receptor cell concentrations to evaluate
cell_concentrations = [0.010 * ureg.millimolar, 0.020 * ureg.millimolar, 0.040 * ureg.millimolar]

# Protocol for 'control' titrations (water-water, buffer-buffer,
# titrations into buffer, etc.)

control_protocol = ITCProtocol(
    'control_protocol',
    sample_prep_method='Plates Quick.setup', # thorough cleaning
    itc_method='ChoderaWaterWater.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
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

sim_protocol = ITCProtocol(
    'SIM protocol',
    sample_prep_method='Plates Quick.setup',
    itc_method='1 Injection SIM.inj',
    analysis_method='Onesite',
    experimental_conditions=dict(target_temperature=25, equilibration_time=300, stir_rate=1000, reference_power=5),
    # TODO: Need to adjust this
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
        10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )

# Protocol for cleaning protocol
cleaning_protocol = ITCProtocol(
    'cleaning protocol',
    sample_prep_method='Plates Clean.setup', # thorough cleaning
    itc_method='ChoderaWaterWater.inj',
    analysis_method='Control',
    experimental_conditions=dict(target_temperature=25, equilibration_time=60, stir_rate=1000, reference_power=5),
    injections=[dict(volume_inj=0.2, duration_inj=0.4, spacing=60, filter_period=0.5)] +
    10 * [dict(volume_inj=3.0, duration_inj=6, spacing=120, filter_period=0.5)],
    )

# Define ITC Experiment.

# use specified protocol by default
itc_experiment_set = ITCExperimentSet(name='Receptor-ligand experiments')
# Add available plates for experiments.
destination_plate = Labware(
    RackLabel='DestinationPlate',
    RackType='ITC Plate')
itc_experiment_set.addDestinationPlate(destination_plate)


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
    name = f'water into water'
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer'
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=buffer_trough,
            cell_source=buffer_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# ligands/receptor
# scale cell concentration to fix necessary syringe concentrations

for cell_concentration in cell_concentrations:
    # Buffer into receptor
    name = f'buffer into receptor'
    experiment = ITCExperiment(
        name=name,
        syringe_source=buffer_trough,
        cell_source=receptor_solution,
        protocol=binding_protocol,
        cell_concentration=cell_concentration,
        buffer_source=buffer_trough,
        cell_volume=cell_volume)
    itc_experiment_set.addExperiment(experiment)

    for ligand, ligand_solution, ligand_ka in zip(ligands, ligand_solutions, ligand_kas):
        # Perform specified number of replicates
        for replicate in range(nreplicates):

            # SIM ligand into buffer with highest concentration of ligand
            name = f'SIM {ligand.name} into buffer (replicate {replicate+1})'
            experiment = ITCExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=buffer_trough,
                protocol=sim_protocol,
                cell_concentration=None,
                buffer_source=buffer_trough,
                cell_volume=cell_volume)
            itc_experiment_set.addExperiment(experiment)

            # SIM ligand into receptor with highest concetration of ligand
            name = f'SIM {ligand.name} into receptor (replicate {replicate+1})'
            experiment = ITCExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=receptor_solution,
                protocol=sim_protocol,
                cell_concentration=cell_concentration,
                buffer_source=buffer_trough,
                cell_volume=cell_volume)
            itc_experiment_set.addExperiment(experiment)

            # Titrate ligand into protein experiment based on guess of affinity
            name = f'titration of {ligand.name} into {receptor.name} (replicate {replicate+1})'
            receptor_experiment = ITCHeuristicExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=receptor_solution,
                protocol=binding_protocol,
                cell_concentration=cell_concentration,
                buffer_source=buffer_trough,
                cell_volume=cell_volume)
            # Optimize the syringe_concentration using heuristic equations and known binding constants
            n_injections = 10 # TODO: Extract this from binding_protocol object
            receptor_experiment.heuristic_syringe(ligand_ka, n_injections, strict=False)

            ligand_experiment = copy.deepcopy(receptor_experiment)
            ligand_experiment.name = f'titration of {ligand.name} into buffer (replicate {replicate+1})'
            ligand_experiment.cell_source = buffer_trough
            ligand_experiment.cell_concentration = None

            # Add experiemnt
            itc_experiment_set.addExperiment(ligand_experiment)
            itc_experiment_set.addExperiment(receptor_experiment)

# Add cleaning experiment.
name = 'final cleaning water titration'
ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=cleaning_protocol, cell_volume=cell_volume)
itc_experiment_set.addExperiment(experiment)

# Water control titrations.
nfinal = 1
for replicate in range(nfinal):
    name = f'water into water'
    experiment = ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume)
    itc_experiment_set.addExperiment(experiment)


# Check that the experiment can be carried out using available solutions and plates.
# Also generate Tecan EVO worklists
import sys
itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)

# Write Tecan EVO pipetting operations for both liquid and air LiHas.
worklist_prefix = 'experiments'
itc_experiment_set.writeTecanWorklist(worklist_prefix)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'experiments.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
