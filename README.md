# itctools

[![Build Status](https://travis-ci.org/choderalab/itctools.svg)](https://travis-ci.org/choderalab/itctools)
[![Build status](https://ci.appveyor.com/api/projects/status/eexd1la2ghs518mc?svg=true)](https://ci.appveyor.com/project/bas-rustenburg/itctools)
[![Coverage Status](https://coveralls.io/repos/choderalab/itctools/badge.png)](https://coveralls.io/r/choderalab/itctools)
[![Code Health](https://landscape.io/github/choderalab/itctools/master/landscape.svg)](https://landscape.io/github/choderalab/itctools/master)

Tools for setting up ITC experiments in an automated fashion using the Tecan EVO and Auto-iTC 200.

## New:

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

You can then install the library using:
```shell
python setup.py install
```

## Requires
-
`anaconda` or `miniconda` with `python` 2.7 or 3.3
Using `conda`, install:
`setuptools`
`numpy`
`openpyxl`
`coverage`
`pip` - to install further dependencies.

Using `pip`, install:
`pint`
`behave`
`coveralls`

## Contains

### `itctools/`

A library for setting up experiments

  - `materials.py`
    - Contains the objects to define chemical compounds, solutions and mixtures
  - `procedures.py`
    - Contains classes that set up experimental procedures, like binding experiments or heat of mixing.
  - `labware.py`
    - Defines chemical containers and their locations on the EVO Deck.
  - `itctools.py`
    - Contains useful functions that don't fall into any other category.

### `examples/`

Some example scripts that use the library.

#### `examples/host_guest`

  - `host_guest.py` prepares a worklist and xlsx file for titrating a host compound into several guests.

#### `examples/mixture_heats`

  -  `mixture_heats.py` prepares a worklist and xlsx file for performing heat of mixing experiments.

### `tests/`

Tests that are used to maintain the integrity of the library.

#### `tests/behave`

Tests that use `behave` to make sure the example scripts are functional. Run `behave` from the root directory of the repository to run these tests. See also the `.behaverc` file.

##### `tests/behave/features`

- The feature files written in `Gherkin` that test the behavior of the example scripts.
- `environment.py` which sets up a temporary directory called `tmp` in your current directory to run tests in. Note that it automatically erases this directory before every run to ensure a clean working directory.

##### `tests/behave/features/steps`
- `scripts.py` holds the step definitions that are run by behave when testing the features defined in the feature files.
