"""
Script for generation of binary mixture ITC titrations.
"""


from simtk.unit import *

from itctools.protocols import ITCProtocol, ITCExperimentSet, ITCExperiment, ITCHeuristicExperiment
from itctools.chemicals import Compound, Solvent, SimpleSolution
from itctools.labware import Labware, PipettingLocation
 


#TODO command line specification of density and name
# Mimic command line input by manually setting variables for now

#Name for the entire set of experiments
set_name = 'water-dmso mixtures'

#Do a cleaning titration at the end of the set
final_cleaning = False

#For every liquid that is to be mixted, define
label1 = 'Water' #identifier
dens1 = 0.9970479 *grams/milliliter #density at standard conditions
mw1= 18.01528*gram / mole    # molecular weight
pur1= 1.0 #purity ( TODO might make this optional)

label2 = 'Dimethyl sulfoxide'
dens2 = 1.092*grams/milliliter
mw2=  78.13 * gram /mole 
pur2= 1.0

#Which liquid to use for controls (Probably water)
#TODO We may want to do the controls with every liquid
control_index = 0


# END of user input


#Define all the liquids used
liquids = list()
liquids.append(PureLiquid(label1, dens1, mw1, purity=pur1))
liquids.append(PureLiquid(label2, dens2, mw2, purity=pur2))
control_liquid = liquids[control_index] #TODO  more than 1 control liquid

# Define the vial holder
source_plate = Labware(RackLabel='SourcePlate', RackType='5x3 Vial Holder')

#Define the location of all liquids in the vial holder
# NOTE Vials must be put in vial holder in the order that they were specified as input
locations = list()
for l,liquid in enumerate(liquids):
    locations.append(PipettingLocation(source_plate.RackLabel, source_plate.RackType, l))
    


# Define Mixing protocols.

# Protocol for 'control' titrations (water-water)
control_protocol = MixingProtocol('control protocol', sample_prep_method='Plates Quick.setup', itc_method='ChoderaWaterWater.inj', analysis_method='Control')
# Protocol for a titration with increasing mole fraction
#TODO Define the mixing protocol at the ITC machine
mixing_protocol = MixingProtocol('mixture protocol',  sample_prep_method='Plates Quick.setup', itc_method='ChoderaHostGuest.inj', analysis_method='Onesite')
# Protocol for cleaning protocol
cleaning_protocol = MixingProtocol('cleaning protocol', sample_prep_method='Plates Clean.setup', itc_method='water5inj.inj', analysis_method='Control')

# Define the experiment set.
mixing_experiment_set = MixingExperimentSet(name=set_name) # use specified protocol by default

# Add available plates for experiments.
mixing_experiment_set.addDestinationPlate(Labware(RackLabel='DestinationPlate', RackType='ITC Plate'))
mixing_experiment_set.addDestinationPlate(Labware(RackLabel='DestinationPlate2', RackType='ITC Plate'))

nreplicates = 1 # number of replicates of each experiment
ncontrols = 1 #initial controls
nfinal = 1 # final (water-water) controls

control_mixture = (components=[control_liquid], molefractions=[1.0], locations=[locations[control_index]], normalize_fractions=False)

# Add cleaning titration

name = 'initial cleaning water titration'
mixing_experiment_set.addExperiment(HeatOfMixingExperiment(name=name, control_mixture, control_mixture, cleaning_protocol))

# Add control titrations.
#TODO Perform control for liquid x into x, for every input liquid?
for replicate in range(ncontrols):
    name = 'water into water %d' % (replicate+1)
    mixing_experiment_set.addExperiment( HeatOfMixingExperiment(name=name, control_mixture, control_mixture, protocol=control_protocol) )

#Define mixing experiments here


#Add cleaning experiment.
if final_cleaning:
    name = 'initial cleaning water titration'
    mixing_experiment_set.addExperiment(HeatOfMixingExperiment(name=name, control_mixture, control_mixture, cleaning_protocol))


# Water control titrations.
# Add control titrations.
#TODO Perform control for liquid x into x, for every input liquid?
for replicate in range(nfinal):
    name = 'water into water %d' % (replicate+1)
    mixing_experiment_set.addExperiment( HeatOfMixingExperiment(name=name, control_mixture, control_mixture, protocol=control_protocol) )

# Check that the experiment can be carried out using available solutions and plates.
#TODO make validation function complete
#itc_experiment_set.validate(print_volumes=True, omit_zeroes=True)


#Allocate experiment resources, destinations, et cetera
mixing_experiment_set.setup_mixing_experiments()


#For convenience, concentrations

# Write Tecan EVO pipetting operations.
#worklist_filename = 'mixing-itc.gwl'
#mixing_experiment_set.writeTecanWorklist(worklist_filename)

# Write Auto iTC-200 experiment spreadsheet.
#excel_filename = 'run-itc.xlsx'
#mixing_experiment_set.writeAutoITCExcel(excel_filename)
