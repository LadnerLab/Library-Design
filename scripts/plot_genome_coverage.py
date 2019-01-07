#!/usr/bin/env python3
from Bio import Entrez           # for efetch
import argparse                  # for parsing command-line arguments
from time import sleep           # For limiting requests
import matplotlib                # for generating charts
matplotlib.use( 'Agg' )          # Don't want to display charts
import matplotlib.pyplot as plt  # For generating charts
import sys                       # for handling errors
import os                        # for directory operations

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

    # parse accession-taxId file




class FileParser:
    def __init__( self, filename = None ):
        self._filename    = filename
        self._parsed_data = None

    def parse( self ):
        pass

class AccessionParser( FileParser ):
    def __init__( self, filename = None ):
        super().__init__( filename )

    def parse( self, filename = None ):
        parse_file = ""
        if filename:
            parse_file = filename
        else:
            parse_file = self._filename

        return self._parse_file( parse_file )

    def _parse_file( self, filename ):
        data_dict = {}
        with open( filename, 'r' ) as open_file:
            for lineno, line in enumerate( open_file ):
                if lineno: # Skip first line
                    split_line = line.strip().split( '\t' )
                    if split_line[ 0 ] not in data_dict:
                        data_dict[ split_line[ 0 ] ] = set()

                    data_dict[ split_line[ 0 ] ].add( split_line[ 1 ] )
        return data_dict
                    
                        
            
        

            
if __name__ == '__main__':
    main()
