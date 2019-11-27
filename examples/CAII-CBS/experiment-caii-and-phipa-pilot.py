#! /usr/bin/env python
"""
Script for generation of input files for CAII:CBS ITC.
"""

from itctools.itctools import ureg, Quantity
from itctools.procedures import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.materials import Solvent, Compound, SimpleSolution
from itctools.labware import Labware, PipettingLocation
import pprint
import copy


# The sample cell volume in microliters
cell_volume = 202.8

# Define solvents.
water = Solvent('water', density=0.9970479 * ureg.gram / ureg.milliliter)
caii_buffer = Solvent('caii_buffer', density=1.014 * ureg.gram / ureg.milliliter) # TODO is our density the same as the HOST-GUEST buffer?
PHIPA_buffer = Solvent('PHIPA_buffer', density=1.014 * ureg.gram / ureg.milliliter) # TODO is our density the same as the HOST-GUEST buffer?

# Define compounds.
caii = Compound('CAII', molecular_weight=29246.0 * (ureg.gram / ureg.mole), purity=1.)
cbs = Compound('CBS', molecular_weight=201.2 * (ureg.gram / ureg.mole), purity=0.97)

PHIPA = Compound('PHIPA', molecular_weight=14117.0 * (ureg.gram / ureg.mole), purity=1.)
cacl = Compound('CaCl2*2H2O', molecular_weight=147.01 * (ureg.gram / ureg.mole), purity=0.97)
hepes = Compound('HEPES', molecular_weight=238.3012 * (ureg.gram / ureg.mole), purity=0.97)

#Ka (association constants) TODO Add this to the compound properties? (maybe a dict with protein as key)
cbs_ka = 1.203e6 / ureg.molar

cacl_ka = 1.0 / (10 * ureg.micromolar)
hepes_ka = 1.0 / (100 * ureg.micromolar)

# Define troughs on the instrument
water_trough = Labware(RackLabel='Water', RackType='Trough 100ml')
caii_buffer_trough = Labware(RackLabel='Buffer', RackType='Trough 100ml')
PHIPA_buffer_trough = Labware(RackLabel='SecondBuffer', RackType='Trough 100ml')

# Define source labware.
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')
protein_source_plate = Labware(RackLabel='ProteinSourcePlate', RackType='ITC Plate')

# Define source solutions in the vial holder
# TODO : Define solutions once prepared with the Quantos
#caii_solution = SimpleSolution(compound=caii, compound_mass=5 * ureg.milligram, solvent=caii_buffer, solvent_mass=0.7 * ureg.gram, location=PipettingLocation(
caii_solution = SimpleSolution(compound=caii, compound_mass=2.5 * ureg.milligram, solvent=caii_buffer, solvent_mass=1.5 * ureg.gram, location=PipettingLocation(
    protein_source_plate.RackLabel,
    protein_source_plate.RackType,
    1)) # Well A1 of ITC plate

# 124 uM @ 1.5 mL
PHIPA_solution = SimpleSolution(compound=PHIPA, compound_mass=2.6 * ureg.milligram, solvent=PHIPA_buffer, solvent_mass=1.5 * ureg.gram, location=PipettingLocation(
    protein_source_plate.RackLabel,
    protein_source_plate.RackType,
    2)) # Well B1 of ITC plate

cbs_solution = SimpleSolution(compound=cbs, compound_mass=4.92 * ureg.milligram, solvent=PHIPA_buffer, solvent_mass=12.954 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    1)) # Well A1 of vial holder

cacl_solution = SimpleSolution(compound=cacl, compound_mass=20.0 * ureg.milligram, solvent=PHIPA_buffer, solvent_mass=10.0 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    2)) # Well A1 of vial holder

hepes_solution = SimpleSolution(compound=hepes, compound_mass=70.0 * ureg.milligram, solvent=caii_buffer, solvent_mass=10.0 * ureg.gram, location=PipettingLocation(
    source_plate.RackLabel,
    source_plate.RackType,
    3)) # Well A1 of vial holder

proteins = ['caii', 'PHIPA']

ligands = dict()
ligands['caii'] = [cbs]
ligands['PHIPA'] = [hepes]

ligand_solutions = dict()
ligand_solutions['caii'] = [cbs_solution]
ligand_solutions['PHIPA'] = [hepes_solution]

ligand_kas = dict()
ligand_kas['caii'] = [cbs_ka]
ligand_kas['PHIPA'] = [hepes_ka]

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
itc_experiment_set = ITCExperimentSet(name='Bovine carbonic anhydrase (CAII) and PHIPA pilot experiments')
# Add available plates for experiments.
destination_plate = Labware(
    RackLabel='DestinationPlate',
    RackType='ITC Plate')
itc_experiment_set.addDestinationPlate(destination_plate)

nreplicates = 1  # number of replicates of each experiment

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

# Add CAII buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=caii_buffer_trough,
            cell_source=caii_buffer_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# buffer into CAII
for replicate in range(1):
    name = 'buffer into CAII %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=caii_buffer_trough,
            cell_source=caii_solution,
            protocol=control_protocol,
            cell_concentration=10 * ureg.micromolar,
            buffer_source=caii_buffer_trough,
            cell_volume=cell_volume))

# ligands/CAII
# scale cell concentration to fix necessary syringe concentrations

protein = 'caii'
for ligand, ligand_solution, ligand_ka in zip(ligands[protein], ligand_solutions[protein], ligand_kas[protein]):

    # We need to store the experiments before adding them to the set
    ligand_protein_experiments = list()
    ligand_buffer_experiments = list()

    # Scaling factors per replicate
    factors = list()

    # Define ligand to protein experiments
    for cell_concentration in [10 * ureg.micromolar, 20 * ureg.micromolar, 40 * ureg.micromolar]:
        for replicate in range(1):
            name = '%s into CAII %d' % (ligand.name, replicate + 1 )
            experiment = ITCHeuristicExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=caii_solution,
                protocol=binding_protocol,
                cell_concentration=cell_concentration,
                buffer_source=caii_buffer_trough,
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
        experiment.cell_source = caii_buffer_trough
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
nfinal = 1
for replicate in range(nfinal):
    name = 'final water into water test %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=water_trough,
            cell_source=water_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

#
# PHIPA
#

# Add CAII buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=PHIPA_buffer_trough,
            cell_source=PHIPA_buffer_trough,
            protocol=control_protocol,
            cell_volume=cell_volume))

# buffer into CAII
for replicate in range(1):
    name = 'buffer into PHIPA %d' % (replicate + 1)
    itc_experiment_set.addExperiment(
        ITCExperiment(
            name=name,
            syringe_source=PHIPA_buffer_trough,
            cell_source=PHIPA_solution,
            protocol=control_protocol,
            cell_concentration=10 * ureg.micromolar,
            buffer_source=PHIPA_buffer_trough,
            cell_volume=cell_volume))

# ligands/CAII
# scale cell concentration to fix necessary syringe concentrations

protein = 'PHIPA'
for ligand, ligand_solution, ligand_ka in zip(ligands[protein], ligand_solutions[protein], ligand_kas[protein]):

    # We need to store the experiments before adding them to the set
    ligand_protein_experiments = list()
    ligand_buffer_experiments = list()

    # Scaling factors per replicate
    factors = list()

    # Define ligand to protein experiments
    for cell_concentration in [10 * ureg.micromolar, 20 * ureg.micromolar, 40 * ureg.micromolar]:
        for replicate in range(1):
            name = '%s into PHIP %d' % (ligand.name, replicate + 1 )
            experiment = ITCHeuristicExperiment(
                name=name,
                syringe_source=ligand_solution,
                cell_source=PHIPA_solution,
                protocol=binding_protocol,
                cell_concentration=cell_concentration,
                buffer_source=PHIPA_buffer_trough,
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
        experiment.cell_source = PHIPA_buffer_trough
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
nfinal = 1
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
with open('experiment.txt', 'w') as outfile:
    itc_experiment_set.validate(print_volumes=True, omit_zeroes=True, human_readable_log=outfile)

# For convenience, concentrations
for protein in proteins:
    for ligand_solution in ligand_solutions[protein]:
        print("%12s %.4f mM" % (ligand_solution.name, ligand_solution.concentration / ureg.millimolar ))
    print("%12s %.4f mM" % ('CAII', caii_solution.concentration  / ureg.millimolar ))
    print("%12s %.4f mM" % ('PHIPA', PHIPA_solution.concentration  / ureg.millimolar ))

# Write Tecan EVO pipetting operations for both liquid and air LiHas.
worklist_prefix = 'experiment'
itc_experiment_set.writeTecanWorklist(worklist_prefix)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'experiment.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
