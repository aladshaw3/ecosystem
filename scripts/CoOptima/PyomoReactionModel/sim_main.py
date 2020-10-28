## Script to run a simulation
from reaction_model_framework import ReactionModel
import os, sys, getopt
import time
import os.path
from os import path
import platform

# Create model and run simulation from given file
def execute_model(file):
    result = ReactionModel()
    result.read_yaml_simfile(file)
    result.run_model()
    result.print_to_file()
    #result.Display()

##Define a help message to display
def help_message():
    print()
    print("sim_main.py script command line options...")
    print()
    print("\t-h           Print this message")
    print("\t-i <folder>  Provide the name of the file (with directory) that contains main simulation control file.")
    print()
    print("Example Usage:")
    print()
    print("\tpython sim_main.py -i ethene/ethene_sim.yml")
    print()
    print("\tNOTE: The control file MUST be formatted using the PyYaml library standard.")
    print("\t\t See https://pyyaml.org/wiki/PyYAMLDocumentation for more details.")
    print()
    print("\tNOTE: Usage of this script REQUIRES the 'ipopt' solver packaged with IDAES.")
    print("\t\t See https://idaes-pse.readthedocs.io/en/stable/advanced_user_guide/advanced_install/index.html#advanced-user-installation for installation details.")
    print()

## Change the input argument formatting if the system is 'Windows'
def update_file_name(file):
    new_file = file
    if platform.system() == 'Windows':
        # Remove the first "." character
        if ".\\" in file:
            new_file = ""
            for item in file.split(".\\")[1].split("\\"):
                if "." not in item:
                    new_file+=item+"/"
                else:
                    new_file+=item
    return new_file

##Define the 'main' function
#
#   argv is the list of arguments pass to the script at the command line
#
#   Accepted arguments include...
#
#       -h         ==>  display help information
#
#       -i dir/    ==>   path and name of the folder than contains other folders of data
#
#       -o dir/    ==>   path and name of the folder to place output into
def main(argv):
    input_file = ""
    #Check for valid arguments
    try:
        opts, args = getopt.getopt(argv,"hi:",["input_file="])
    except:
        help_message()
        sys.exit(2)
    #Iterate through arguments and save some information
    for opt, arg in opts:
        if opt == '-h':
            help_message()
            sys.exit()
        if opt in ("-i", "--ifile"):
            input_file = arg

    #If we made it to this point, then no errors in input
    #Check to see if the input_folder does exist
    if path.exists(input_file) == False:
        print("Error! Given argument is not a valid file or location to file!")
        help_message()
        sys.exit()

    start = time.time()
    input_file = update_file_name(input_file)
    execute_model(input_file)
    end = time.time()
    elapse_min = (end-start)/60
    print("\nCOMPLETED!!!")
    print("\tRuntime (min): " + str(elapse_min) + "\n")


##Directs python to call the main function
if __name__ == "__main__":
   main(sys.argv[1:])
