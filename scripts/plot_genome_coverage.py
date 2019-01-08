#!/usr/bin/env python3
from Bio import Entrez           # for efetch
import argparse                  # for parsing command-line arguments
import time                      # for tracking request frequency, sleeping
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

class RecordWriter:
    def __init__( self, suffix = "", prefix = "", work_dir = "" ):
        self._suffix   = suffix
        self._prefix   = prefix
        self._work_dir = work_dir 

        if not work_dir:
            self._work_dir = os.getcwd()

    def write_file( self, filename, new_record, append = False ):
        open_mode = 'w'
        if append:
            open_mode = 'a'

        start_dir = os.getcwd()

        os.chdir( self._work_dir )
        openfile_name = "%s%s%s" % ( self._prefix, filename, self_suffix )

        with open( openfile_name, open_mode ) as open_file:
            open_file.write( new_record.strip() )

        os.chdir( start_dir )
        
class EntrezController:
    def __init__( self, email = None, database = None, api_key = None,
                  rettype = None, retmode = None
                ):
        self._conn = EntrezConnection( email = email,
                                       api_key = api_key
                                     )
        self._database = database
        self._rettype  = rettype
        self._retmode  = retmode

    def set_database( self, new_db ):
        self._database = database

    def set_rettype( self, new_type ):
        self._rettype = new_type

    def set_retmode( self, new_mode ):
        self._retmode = new_mode

    def set_connection( self, new_connection ):
        self._conn = new_connection

    def get_record( self, id ):
        return self._conn.query( Entrez.efetch, db = self._database,
                                 rettype = self._rettype, retmode = self._retmode,
                                 id = id
                               ).read()
                                 
class EntrezConnection:
    API_KEY_REQS_PER_SEC = 10
    EMAIL_REQS_PER_SEC   = 3

    def __init__( self, email = None, api_key = None ):
        self._email               = email
        self._sleep_time          = 1
        self._api_key             = api_key

        Entrez.email              = email
        Entrez.api_key            = api_key 

        self._requests_this_second = 0
        self._time                 = time.time()

        if email:
            self._max_reqs_per_second = EntrezConnection.EMAIL_REQS_PER_SEC
        elif api_key:
            self._max_reqs_per_second = EntrezConnection.API_KEY_REQS_PER_SEC

    class EmailNotSuppliedException( Exception ):
        def __str__( self ):
            return "Email address was not supplied for use with Entrez"

    def set_email( self, new_email ):
        self._email = new_email
        Entrez.email = new_email

    def set_key( self, new_key ):
        self._api_key  = new_key
        Entrez.api_key = new_key

    def set_max_requests_per_second( self, new_reqs ):
        self._max_reqs_per_second = new_reqs

    def query( self, function, **kwargs ):
        if not self._email:
            raise EntrezConnection.EmailNotSuppliedException()

        now = time.time()
        if now - self._time <= 1.0:
            self._requests_this_second += 1

            if self._requests_this_second > \
               self._max_reqs_per_second:
                time.sleep( self._sleep_time )
        else:
            self._requests_this_second = 1

        self._time = now
        return function( **kwargs )

class MissingArgumentException( Exception ):
    def __init__( self, str_reason ):
        self._reason = str_reason
        
    def __str__( self ):
        return "%s was either not provided from the command line, " \
               "or was provided with no argument" % self._reason

if __name__ == '__main__':
    main()
