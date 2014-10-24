itctools
=========
[![Build Status](https://travis-ci.org/choderalab/itctools.svg)](https://travis-ci.org/choderalab/itctools)

Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200.

New:
---
Run `bash install.sh` to set up a new python environment using miniconda at ${HOME}/miniconda. It should take care of all dependencies. 

To access the new python enviroment you can type 
```shell
source ${HOME}/miniconda/bin/activate 2.7
```
or extend your `${PATH}` variable with `${HOME}/miniconda/bin`:
```shell
export PATH=${PATH}:${HOME}/miniconda/bin
```
To use the scripts within an existing environment, please look at the dependencies and install them accordingly. You can use the requirements files from the installer.

Requires
--------
`anaconda` or `miniconda` with `python` 2.7 or 3.3  
`simtk` (which is part of `openmm`.)  
`numpy`  
`openpyxl`  
`pip` - to install:  
`behave`

Contains
--------
`itctools/`

A library for setting up experiments

  - `materials.py`
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
    
`features/`

Test the integrity of the library. Run `behave` from top directory in order to make sure the scripts are functional.
