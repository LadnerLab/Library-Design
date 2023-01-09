#!/usr/bin/env python

# Add or remove modules here
import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules
import glob, os, shutil

from collections import defaultdict


# the main function contains the bulk of the script
def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # These are optional arguments
    parser.add_argument("-f", "--filler", default="XXXXXX", help="Sequence to use as filler for designs.")
    parser.add_argument("-s", "--spacer", default="XXXXXX", help="Sequence to use as spacer for designs.")
    
    # These are arguments that the user is required to provide
    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-i", "--input", help="Base sequences for design. Tab-delimited file with ... .", required=True)

    args = parser.parse_args()
    
    # Add code for script here
    
    
#----------------------End of main()

# Put any other function definitions here


###------------------------------------->>>>    

# This allows the script to be imported as a module to other scripts, without running the main() code
# Therefore, you could use the other defined functions in another script, if you wanted
if __name__ == "__main__":
    main()

