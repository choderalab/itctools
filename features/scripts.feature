Feature: Heat of mixing experiments
  In order to perform heat of mixing experiments using the AutoITC200,
  a tecan EVO worklist (.gwl) file is generated to prepare the plate
  To run the ITC experiment itself, a worksheet (.xls) file needs to be generated
  This file tells the AutoITC200 machine which wells to use as cell and syringe solutions
  It sets up wells with various mole fractions from 0 to 1
  The current version uses DMSO and water as examples, but is known not to produce good results
  Reasons for this are the large heats generated because of the large mismatch in concentrations

  Scenario: Running the script
    Given a directory called mixture-itc
     When the python script mixture_heats.py is called
     Then a file called mixing-itc.gwl is created
      And mixing-itc.gwl is formatted as a .gwl file
      And a file called mixing-itc.xlsx is created
      And mixing-itc.xlsx is formatted as a .xslx file


Feature: Host-Guest binding ITC experiment
    Generate a Tecan EVO worklist (.gwl) file to prepare a plate for Host-Guest ITC experiments
    Also generates a worksheet (.xlsx) file containing instructions for the AutoITC200

    Scenario: Running the script
      Given a directory called host-guest-itc
       When the python script host_guest.py is called
       Then a file called host-guest-itc.gwl is created
        And host-guest-itc.gwl is formatted as a .gwl file
        And a file called host-guest-itc.xlsx is created
        And host-guest-itc.xlsx is formatted as a .xslx file