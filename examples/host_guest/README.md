# Host Guest Example

This script produces the input file for measuring the binding affinity of a series of guest compounds to the host cucurbit[7]uril. It uses known estimates of the binding affinity in a heuristic equation (`reference missing`) to predict optimal concentrations to perform the experiment. 

## output

- `host-guest.gwl`, a Tecan EVO worklist that produces buffer solutions for the host and guest compounds at the desired concentrations in a 96 well plate.
- `host-guest.xlsx`, a spreadsheet for the Auto-iTC200 that defines the cell and syringe sources for each titration series. 
