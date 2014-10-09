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
