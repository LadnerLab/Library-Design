#!/usr/bin/env python3
from Bio import Entrez           # for efetch
import argparse                  # for parsing command-line arguments
import time                      # for tracking request frequency, sleeping
import matplotlib                # for generating charts
matplotlib.use( 'Agg' )          # Don't want to display charts
import matplotlib.pyplot as plt  # For generating charts
import sys                       # for handling errors
import os                        # for directory operations
import subprocess                # for running queries

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
    arg_parser.add_argument( '--dont_download', help = "Include this flag if the sequences will not need to be downloaded from "
                                                       "Entrez. If this flag is included, the files will be searched for locally.",
                             action = 'store_true'
                           )
    arg_parser.add_argument( '-t', '--num_threads', help = "Number of threads to use for blast+ queries",
                             default = 1, type = int
                           )
    arg_parser.add_argument( '-l', '--library', help = "Filename in which a designed library is stored, must be in FASTA format. "
                                                       "Note that a blast db will be created from this fasta to improve blasting performance"
                           )

    # parse args
    args = arg_parser.parse_args()

    BLAST_OUT_DIR = 'blast_results'

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
    to_download = not args.dont_download

    # Everyone shares a connection so request limit is not overriden
    controllers.call( EntrezController.set_connection, connection )

    # download genome information, if necessary. 
    if to_download:
        for record in accession_data.as_list():
            filename = record.get_id()
            accession = ','.join( record.get_accession_num() )
            new_records = controllers.call( EntrezController.get_record,
                                            accession
                                          )

            for index, record in enumerate( new_records ):
                writer = writers.as_list()[ index ]
                writer.write_file( filename, record, append = True)

    blaster = SubprocessRunner()

    create_blast_db( blaster, args.library )

    # create a directory to hold the blast outputs 
    if not os.path.exists( BLAST_OUT_DIR ):
        os.mkdir( BLAST_OUT_DIR )

    
    # perform the blast analyses on each protein sequence
    for record in accession_data.as_list():

        # create the command to be used for blasting
        blast_command = make_blast_command( record, args.library,
                                            protein_writer.get_suffix(),
                                            protein_writer.get_work_dir(),
                                            args.num_threads
                                          )


        # invoke the command with the correct input file
        blaster.invoke( blast_command )

    # combine the blast outputs to a single output


    # create a plot, one for each entry in the protein dir

def create_blast_db( runner, ref_file ):
    command = [ 'makeblastdb', '-in %s' % ref_file,
                '-input_type fasta',
                '-dbtype prot'
              ]

    runner.invoke( command )

def make_blast_command( record_data, db_name, file_suffix, work_dir, num_threads ):
    blast_command = [ 'blastp',
                      '-query %s' % ( '%s/%s%s' % ( work_dir,
                                                    record_data.get_id(),
                                                    file_suffix
                                                  )
                                    ),
                      '-db %s' % ( db_name ),
                      '-outfmt 5', # 5 for XML format
                      '-num_threads %d' % num_threads
                      ]

    return blast_command

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
            
class FileParser:
    def __init__( self, filename = None ):
        self._filename    = filename
        self._parsed_data = None

    def parse( self ):
        pass

class SubprocessRunner():
    def __init__( self ):
        self._runner       = subprocess.check_output
        self._last_command = None
        self._error        = False

    def invoke( self, command_arg_list ):
        self._last_command = command_arg_list

        try:
            value = self._runner( command_arg_list ).decode()
            # If we reach this point, runner returned exit status of zero
            self._error = False
            return value

        except:
            self._error = True
            return None

    def get_last_command( self ):
        return self._last_command

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

    def get_suffix( self ):
        return self._suffix
    def get_work_dir( self ):
        return self._work_dir
        
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
