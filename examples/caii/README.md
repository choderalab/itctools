# Carbonic anhydrase II : caboxybenzenesulfonamide (CAII:CBS)

This script produces the input file for measuring the binding affinity of CBS to CAII

## Reference

* Sarah E. Boyce, Joel Tellinghuisen, and John D. Chodera.
  Avoiding accuracy-limiting pitfalls in the study of protein-ligand interactions with isothermal titration calorimetry.
  http://dx.doi.org/10.1101/023796

## Manifest

* `experiment.py` - script to generate experiment worklists and Excel spreadsheet
* `caii.gwl` - a Tecan EVO worklist that produces buffer solutions for the CAII and CBS at the desired concentrations in a 96 well plate
* `caii.xlsx` - a spreadsheet for the Auto-iTC200 that defines the cell and syringe sources for each titration series
