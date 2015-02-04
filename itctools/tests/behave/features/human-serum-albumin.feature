Feature: Human Serum Albumin binding ITC experiment
    Generate a Tecan EVO worklist (.gwl) file to prepare a plate for Human Serum Albumin ITC experiments with aspirin and naproxen
    Also generates a worksheet (.xlsx) file containing instructions for the AutoITC200

    Scenario: Running the script
      Given that the script is in the directory "examples/human-serum-albumin"
        And the module "itctools" is installed 
        And the working directory is "human-serum-albumin"
       When the script "hsa.py" is called successfully from the working directory
       Then a file called "hsa.gwl" is created
        And "hsa.gwl" is not an empty file
        And a file called "hsa.xlsx" is created
        And "hsa.xlsx" is not an empty file
