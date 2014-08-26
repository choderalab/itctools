itctools
=========
Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200.

Requires
--------
python 2.7 or above

simtk

itertools

nose

Contains
--------
`itctools/`

A library for setting up experiments

  - `chemicals.py`
    - Contains the objects to define chemical compounds, solutions and mixtures
  - `procedures.py`
    - Contains classes that set up experimental procedures, like binding experiments or heat of mixing.
  - `labware.py`
    - Defines chemical containers and their locations on the EVO Deck.
  - `itctools.py`
    - Contains useful functions that don't fall into any other category. 

`scripts/`

Some example scripts that use the library.

  - `host_guest.py`
    - Prepares a worklist and xlsx file for titrating a host compound into several guests.
  - `mixture_heats.py`
    - Prepare a worklist and xlsx file for performing heat of mixing experiments.
    
`tests.py`

Test the integrity of the library. Run `nosetests` in order to make sure the script is functional.
