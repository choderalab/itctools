"""
Script for generation of binary mixture ITC titrations.
"""


from simtk.unit import *

from itctools.protocols import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.chemicals import Solvent, SimpleSolution
from itctools.labware import Labware, PipettingLocation
 


#TODO command line specification of density and name
# Mimic command line input by manually setting variables for now
label1 = 'Water'
dens1 = 0.9970479 *grams/milliliter

label2 = 'Dimethyl sulfoxide'
dens2 = 1.092*grams/milliliter

# Define mixture components
comp1 = Solvent(label1, density=dens1)
comp2 = Solvent(label2, density=dens2)

# Not using troughs. #
#Definition of an example trough
#comp1_trough = Labware(RackLabel=label1, RackType='Trough 100ml')


# Define source labware
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')

# Define source solutions on the deck.one

host_solution = SimpleSolution(compound=host, compound_mass=16.76*milligrams, solvent=buffer, solvent_mass=10.2628*grams, location=PipettingLocation(source_plate.RackLabel, source_plate.RackType, 1))
guest_solutions = list()

# Define ITC protocolss.

# Protocol for 'control' titrations (water-water)
control_protocol = ITCProtocol('control protocol', sample_prep_method='Plates Quick.setup', itc_method='ChoderaWaterWater.inj', analysis_method='Control')
# Protocol for 1:1 binding analyis
mixing_protocol = ITCProtocol('mixture protocol',  sample_prep_method='Plates Quick.setup', itc_method='ChoderaHostGuest.inj', analysis_method='Onesite')
# Protocol for cleaning protocol
cleaning_protocol = ITCProtocol('cleaning protocol', sample_prep_method='Plates Clean.setup', itc_method='water5inj.inj', analysis_method='Control')

# Define ITC Experiment.

itc_experiment_set = ITCExperimentSet(name='SAMPL4-CB7 host-guest experiments') # use specified protocol by default
# Add available plates for experiments.
itc_experiment_set.addDestinationPlate(Labware(RackLabel='DestinationPlate', RackType='ITC Plate'))
itc_experiment_set.addDestinationPlate(Labware(RackLabel='DestinationPlate2', RackType='ITC Plate'))

nreplicates = 1 # number of replicates of each experiment

# Add cleaning experiment.
name = 'initial cleaning water titration'
itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=cleaning_protocol) )

# Add water control titrations.
for replicate in range(1):
    name = 'water into water %d' % (replicate+1)
    itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=control_protocol) )

# Add buffer control titrations.
for replicate in range(1):
    name = 'buffer into buffer %d' % (replicate+1)
    itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=buffer_trough, cell_source=buffer_trough, protocol=control_protocol) )

# Host into buffer.
for replicate in range(1):
    name = 'host into buffer %d' % (replicate+1)
    itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=host_solution, cell_source=buffer_trough, protocol=binding_protocol) )

# Host/guests.
# scale cell concentration to fix necessary syringe concentrations
cell_scaling = 1.
for guest_index in range(nguests):

    #We need to store the experiments before adding them to the set
    host_guest_experiments = list()
    buff_guest_experiments = list()

    #Scaling factors per replicate
    factors = list()

    # Define host into guest experiments.
    for replicate in range(1):
        name = 'host into %s' % guests[guest_index].name
        experiment = ITCHeuristicExperiment(name=name, syringe_source=host_solution, cell_source=guest_solutions[guest_index], protocol=binding_protocol, cell_concentration=0.2*millimolar*cell_scaling, buffer_source=buffer_trough)
        #optimize the syringe_concentration using heuristic equations and known binding constants
        #TODO extract m, v and V0 from protocol somehow?
        experiment.heuristic_syringe(guest_compound_Ka[guest_index], 10, 3. * microliters, 202.8 * microliters)
        #rescale if syringe > stock. Store factor.
        factors.append(experiment.rescale())
        host_guest_experiments.append(experiment)

    # Define buffer into guest experiments.
    for replicate in range(1):
        name = 'buffer into %s' % guests[guest_index].name
        experiment = ITCHeuristicExperiment(name=name, syringe_source=buffer_trough, cell_source=guest_solutions[guest_index], protocol=blank_protocol, cell_concentration=0.2*millimolar, buffer_source=buffer_trough)
        #rescale to match host into guest experiment concentrations.
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
    name = 'final water into water test %d' % (replicate+1)
    itc_experiment_set.addExperiment( ITCExperiment(name=name, syringe_source=water_trough, cell_source=water_trough, protocol=control_protocol) )

# Check that the experiment can be carried out using available solutions and plates.

itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)

#For convenience, concentrations
for g,guest in enumerate(guest_solutions, start=1):
    print "guest%02d"% g, guest.concentration.in_units_of(millimolar)

print "host", host_solution.concentration.in_units_of(millimolar)


# Write Tecan EVO pipetting operations.
worklist_filename = 'setup-itc.gwl'
itc_experiment_set.writeTecanWorklist(worklist_filename)

# Write Auto iTC-200 experiment spreadsheet.
excel_filename = 'run-itc.xlsx'
itc_experiment_set.writeAutoITCExcel(excel_filename)
