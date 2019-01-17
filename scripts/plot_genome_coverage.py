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
    arg_parser.add_argument( '-e', '--email', help = "Email address to use for creating an Entrez connection, "
                                                     "either this argument, or api_key must be supplied so that "
                                                     "a connection can be made. With just an email, only 3 requests "
                                                     "will be made per second."
                           )
    arg_parser.add_argument( '--api_key', help = "Api key to use for an Entrez connection. "
                                                 "With a valid api key, 10 requests will be made "
                                                 "per second."
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

    connection        = EntrezConnection( email = args.email, api_key = args.api_key )
    genome_controller     = EntrezController( database = args.database, rettype = 'fasta', retmode = 'text' )
    nucleotide_controller = EntrezController( database = args.database, rettype = 'fasta_cds_na', retmode = 'text' )
    protein_controller    = EntrezController( database = args.database, rettype = 'fasta_cds_aa', retmode = 'text' )

    genome_writer         = RecordWriter( suffix = ".fasta", work_dir = "genome_info" )
    nucleotide_writer     = RecordWriter( suffix = ".fna",   work_dir = "nucleotide_info" )
    protein_writer        = RecordWriter( suffix = ".faa",   work_dir = "protein_info" )

    controllers = Composite( genome_controller,
                             nucleotide_controller,
                             protein_controller
    )
    writers    = Composite( genome_writer,
                            nucleotide_writer,
                            protein_writer
                          )

    # Everyone shares a connection so request limit is not overriden
    controllers.call( EntrezController.set_connection, connection )

    # download genome information, if necessary. 
   for record in accession_data.as_list():
       filename = record.get_id()
       accession = ','.join( record.get_accession_num() )
       new_records = controllers.call( EntrezController.get_record,
                                       accession
                                     )

       for index, record in enumerate( new_records ):
           writer = writers.as_list()[ index ]

           writer.write_file( filename, record, append = True)
        
class Composite:
    def __init__( self, *args ):
        self._items = list()
        for current in args:
            self._items.append( current )
            
    def call( self, function, *args, **kwargs ):
        results = list()
        for item in self._items:
            results.append( function( item, *args, **kwargs ) )
        return results

    def call_positional( self, function, arg, **kwargs ):
        results = list()

        assert len( arg ) == len( self._items )

        for index, item in enumerate( self._items ):
            results.append( function( item, args[ index ], **kwargs ) )
        return results

    def as_list( self ):
        return self._items
            

class RecordRetriever:
    def __init__( self,  retriever = None ):
        self._retriever = retriever

    def retrieve( self, record_id ):
        pass

    def set_retriever( self, new_retriever ):
        self._retriever = new_retriever
  
class LocalRecordRetriever( RecordRetriever ):
    def __init__( self, location = None ):
        self._location = location

    def set_location( self, new_loc ):
        self._location = new_loc

    def set_retrieve_method( self, new_method ):
        self._get_record = new_method

    class RecordNotFoundException( Exception ):
        def __init__( self, record ):
            self._not_found = record
            
        def __str__( self ):
            return "%s was not found" % self._not_found

    def retrieve( self, record_id ):
        found = False
        for item in os.listdir( self._location ):
            no_ext = self._remove_file_extension( item )
            if no_ext == record_id:
                found = True

                return_val = [ self._get_record( self._location + '/' + item ) ]

        if not found:
            raise LocalRecordRetriever.RecordNotFoundException( record_id )

        return return_val

    def _get_record( self, filename ):
        record_path = filename
        out_lines = list()

        with open( record_path, 'r' ) as open_file:
            for line in open_file:
                out_lines.append( line.strip() )
        return out_lines
                

    def _remove_file_extension( self, item ):
        return item.strip().split( '.' )[ 0 ]

        
class RemoteRecordRetriever( RecordRetriever ):
    def __init__( self,  retriever = None ):
        super.__init__( retriever )
        
    def retrieve( self, record_id ):
        return self._retriever.get_record( record_id )
        
class DualRecordRetriever:
    def __init__( self, local_retr = None, remote_retr = None ):
        super.__init( retriever )

        self._local_retr  = local_retr
        self._remote_retr = remote_retr

    def retrieve( self, record_id ):
        try:
             yield self._local_retr.retrieve( record_id )
        except LocalRecordRetriever.RecordNotFoundException:
            yield self._remote_retr.retrieve( record_id )

    def set_local_retriever( self, new_retr ):
        self._local_retr = new_retr
        
    def set_remote_retriever( self, new_retr ):
        self._remote_retr = new_retr

    def set_retrievers( self, local_retr = None, remote_retr = None ):
        self.set_local_retriever( local_retr )
        self.set_remote_retriever( remote_retr )
            
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
                        line_data = self._parse_line( line, data_dict )
                    except ValueError:
                        raise AccessionParser.FormatException( lineno, line )
                        
        return data_dict

    def _parse_line( self, string_line, data_dict ):
        split_line = string_line.strip().split( '\t' )
        self._validate_line_format( split_line )

        if split_line[ 0 ] not in data_dict:
            data_dict[ split_line[ 0 ] ] = set()
            
        data_dict[ split_line[ 0 ] ].add( split_line[ 1 ] )

    def _validate_line_format( self, split_line ):
        if len( split_line ) != 2 \
           or len( split_line[ 0 ] ) <= 0 \
              or len( split_line[ 1 ] ) <= 0:
              raise ValueError()
        
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
    def as_list( self ):
        return self._data

class AccessionData:
    def __init__( self, tax_id, accession_num ):
        self._tax_id        = tax_id
        self._accession_num = accession_num 

    def get_id( self ):
        return self._tax_id
    def get_accession_num( self ):
        return self._accession_num

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
        self._create_mode = 0o755

        if not work_dir:
            self._work_dir = os.getcwd()

    def _create_dir( self, dir_name, create_mode ):
        if not os.path.exists( dir_name ):
            os.mkdir( dir_name, create_mode )

    def write_file( self, filename, new_record, append = False ):
        open_mode = 'w+'

        # Only creates the directory if it does not exist
        self._create_dir( self._work_dir, self._create_mode )

        if append:
            open_mode = 'a+'

        start_dir = os.getcwd()

        os.chdir( self._work_dir )
        openfile_name = "%s%s%s" % ( self._prefix, filename, self._suffix )

        with open( openfile_name, open_mode ) as open_file:
            open_file.write( new_record.strip() + '\n' )

        os.chdir( start_dir )

    def set_work_dir( self, new_dir ):
        self._work_dir = new_dir

    def set_suffix( self, new_suffix ):
        self._suffix = new_suffix
        
    def set_prefix( self, new_prefix ):
        self._prefix = new_prefix

class EntrezController:
    def __init__( self, email = None, database = None, api_key = None,
                  rettype = None, retmode = None, create_conn = True
                ):

        if create_conn:
            self._conn = EntrezConnection( email = email,
                                           api_key = api_key
                                         )
        else:
            self._conn = None
        
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

    def get_connection( self ):
        return self._conn

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
        self._max_reqs_per_second = 0

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

        Entrez.email = self._email 
        
        now = time.time()
        if now - self._time <= 1.0:
            self._requests_this_second += 1

            if self._requests_this_second > \
               self._max_reqs_per_second:
                time.sleep( self._sleep_time )
        else:
            self._requests_this_second = 1

        self._time = now
        return_val = None

        added = False

        while not added:
            try:
                return_val = function( **kwargs )
                added = True
            except Exception:
                time.sleep( self._sleep_time )
                self._requests_this_second = 1
            
        return return_val

class MissingArgumentException( Exception ):
    def __init__( self, str_reason ):
        self._reason = str_reason
        
    def __str__( self ):
        return "%s was either not provided from the command line, " \
               "or was provided with no argument" % self._reason

if __name__ == '__main__':
    main()
