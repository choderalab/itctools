Feature: Host-Guest binding ITC experiment
    Generate a Tecan EVO worklist (.gwl) file to prepare a plate for Host-Guest ITC experiments
    Also generates a worksheet (.xlsx) file containing instructions for the AutoITC200

    Scenario: Running the script
      Given that the script is in the directory "examples/host_guest"
        And the module "itctools" is installed 
        And the working directory is "host-guest-itc"
       When the script "host_guest.py" is called successfully from the working directory
       Then a file called "host-guest-itc.gwl" is created
        And "host-guest-itc.gwl" is not an empty file
        And a file called "host-guest-itc.xlsx" is created
        And "host-guest-itc.xlsx" is not an empty file
