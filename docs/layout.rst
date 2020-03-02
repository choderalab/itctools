Library Contents
================

itctools
___________

A library for setting up experiments

  - `materials.py`
    - Contains the objects to define chemical compounds, solutions and mixtures
  - `procedures.py`
    - Contains classes that set up experimental procedures, like binding experiments or heat of mixing.
  - `labware.py`
    - Defines chemical containers and their locations on the EVO Deck.
  - `itctools.py`
    - Contains useful functions that don't fall into any other category.

examples
__________

Some example scripts that use the library.

`examples/host_guest`

  - `host_guest.py` prepares a worklist and xlsx file for titrating a host compound into several guests.

`examples/mixture_heats`

  -  `mixture_heats.py` prepares a worklist and xlsx file for performing heat of mixing experiments.


tests
_______

Tests that are used to maintain the integrity of the library.

 - `tests/behave`
    Tests that use `behave` to make sure the example scripts are functional. Run `behave` from the root directory of the repository to run these tests. See also the `.behaverc` file.

 - `tests/behave/features`
    The feature files written in `Gherkin` that test the behavior of the example scripts.

  - `environment.py` which sets up a temporary directory called `tmp` in your current directory to run tests in. Note that it automatically erases this directory before every run to ensure a clean working directory.

 - `tests/behave/features/steps`
  - `scripts.py` holds the step definitions that are run by behave when testing the features defined in the feature files.
