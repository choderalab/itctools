itctools
=========
Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200.

Requires
--------
python 2.7 or above

simtk

itertools

Contains
--------
itctools/

A library for setting up experiments

  - chemicals.py
    - Contains the objects to define chemical compounds, solutions and mixtures
  - procedures.py
    - Contains classes that set up experimental procedures, like binding experiments or heat of mixing.
  - labware.py
    - Defines chemical containers and their locations on the EVO Deck.

scripts/

Some example scripts that use the library.
