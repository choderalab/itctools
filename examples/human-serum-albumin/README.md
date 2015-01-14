# Human Serum Albumin Example

This script produces the input file for measuring the binding affinity of naproxen and aspirin to human serum albumin (HSA). It uses known estimates of the binding affinity in a heuristic equation (`reference missing`) to predict optimal concentrations to perform the experiment. 

## output

- `hsa.gwl`, a Tecan EVO worklist that produces buffer solutions for the HSA and drugs at the desired concentrations in a 96 well plate.
- `hsa.xlsx`, a spreadsheet for the Auto-iTC200 that defines the cell and syringe sources for each titration series. 
