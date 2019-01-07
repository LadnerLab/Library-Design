#!/usr/bin/env python3
from Bio import Entrez           # for efetch
import argparse                  # for parsing command-line arguments
from time import sleep           # For limiting requests
import matplotlib                # for generating charts
matplotlib.use( 'Agg' )          # Don't want to display charts
import matplotlib.pyplot as plt  # For generating charts

def main():

    # define command-line arguments
    arg_parser = argparse.ArgumentParser( description = "Download genomes given a list of tax ids and "
                                                        "accession numbers, plot BLASTp scores "
                                                        "for each gene."

                                        )
    arg_parser.add_argument( '-a', '--accession', help = "File containing mappings for TaxId's -> "
                                                         "Accession numbers."
                           )

    # parse args
    args = arg_parser.parse_args()


    # verify command-line arguments



if __name__ == '__main__':
    main()
