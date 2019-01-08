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
    arg_parser.add_argument( '-d', '--database', help = "Database from which to retrieve records. "
                                                        "Must be a valid Entrez database name.",
                             default = 'nuccore'
                           )

    # parse args
    args = arg_parser.parse_args()

    # verify command-line arguments
    try:
        verify_args( args )
    except Exception as e:
        ErrorHandler.handle( e, exit = True )
        
    # parse accession-taxId file
    try:
        accession_data = AccessionParser().parse( args.accession )
        accession_data = AccessionDataFactory.from_dict( accession_data )
    except Exception as e:
        ErrorHandler.handle( e, exit = True )

    # download genome information 

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
                    try:
                        line_data = self._parse_line( line )
                    except Exception:
                        # raise AccessionParser.FormatException( lineno, line )
                        pass
                        
        return data_dict

    def _parse_line( self, string_line ):
        split_line = line.strip().split( '\t' )
        self._validate_line_format( split_line )

        if split_line[ 0 ] not in data_dict:
            data_dict[ split_line[ 0 ] ] = set()
            
        data_dict[ split_line[ 0 ] ].add( split_line[ 1 ] )

    def _validate_line_format( self, split_line ):
        if len( split_line ) != 2 \
           or len( split_line[ 0 ] ) <= 0 \
              or len( split_line[ 1 ] ) <= 0:
              raise Exception()
        
    class FormatException( Exception ):
        def __init__( self, line_number, line ):
            self._bad_data = ( line, line_number )
        def __str__( self ):
            return ( "Incorrectly formatted AccessionData: "
                    "%s on line: %d" % ( self._bad_data )
                   )
                    
class Handler:
    def handle( to_handle ):
        pass

class ErrorHandler( Handler ):
    def handle( exception, exit = False ):
        print( "An error has occurred: %s." % str( exception ) )
        if exit:
            print( "Program will exit, "
                   "please fix the above error and restart script."
                 )
            sys.exit( 1 )

class AccessionDataFactory:
    def from_dict( dict_to_get ):
        collection = AccessionDataCollection()
        for tax_id, accession_num in dict_to_get.items():
            collection.add( AccessionData( tax_id, accession_num ) )
        return collection

class AccessionDataCollection:
    def __init__( self ):
        self._data = list()
    def add( self, new_data ):
        self._data.append( new_data )

class AccessionData:
    def __init__( self, tax_id, accession_num ):
        self._tax_id        = tax_id
        self._accession_num = accession_num 

def verify_args( args ):
    if not args.accession:
        raise MissingArgumentException( "TaxID/Accession Map" )
    if not args.database:
        raise MissingArgumentException( "database" )

class MissingArgumentException( Exception ):
    def __init__( self, str_reason ):
        self._reason = str_reason
        
    def __str__( self ):
        return "%s was either not provided from the command line, " \
               "or was provided with no argument" % self._reason

if __name__ == '__main__':
    main()
