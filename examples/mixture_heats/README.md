# Mixture heats example

This example produces input files for a titration series that mixes water and DMSO at various mole fractions to measure heats of mixing. The experiment requires input files for a Tecan EVO and a Microcal Auto-iTC200. Unfortunately, this particular experiment does not produce useful results, but the scripts is kept here for reference and to test the library.

## Output 

- `mixing-itc.gwl`, a Tecan EVO worklist that prepares a 96 well plate with the necessary cell and syringe solutions to perform the experiment.
- `mixing-itc.xlsx`, the corresponding Auto-iTC200 spreadsheet that defines the cell and syringe sources for the titrations that are to be carried out. 
