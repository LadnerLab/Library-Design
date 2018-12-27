#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import subprocess                     # for getting jobstats data
from abc import ABC, abstractmethod   # for abstract base class
import datetime                       # for getting the date and verifying date format
from enum import Enum                 # for error handling
import sys                            # for exiting on abnormal execution
import protein_oligo_library as oligo # for getting kmer-counts


def main():
    # create parser object
    arg_parser = argparse.ArgumentParser( description = "Turn a list of jobid_arrayid|Filename entries into a "
                                                        "list of items containing parseable runtime stats and kmer counts."
                                        )

    # add program options and arguments
    arg_parser.add_argument( '-j', '--job_data', type = str,
                             help = "File containing entries that take the form: "
                                    "filename|'$SLURM_ARRAY_JOB_ID'_'$SLURM_ARRAY_TASK_ID' "
                                    "entries, one per line."
                           )

    arg_parser.add_argument( '-f', '--fasta_dir', type = str,
                             help = "Directory containing fasta files, "
                                    "each filename found in the job_data argument "
                                    "should be found in this directory."
                                    "Note that if this argument is not provided, "
                                    "a kmer_count file will not be produced."
                           )

    arg_parser.add_argument( '-o', '--job_output', type = str, default = 'job_output.txt',
                             help = "Name of file to write output job data to, each entry "
                                    "will take the form: job_id|$data, where data is the jobstats "
                                    "entry for job_id."
                           )

    arg_parser.add_argument( '-k', '--kmer_output', type = str, default = "kmer_count.txt",
                             help = "Name of file to write kmer counts to, "
                                    "one entry will be written per file found in "
                                    "fasta dir, each entry takes the form: "
                                    "filename|num_kmers"
                           )

    arg_parser.add_argument( '-s', '--kmer_size', type = int, default = 9,
                             help = "Size of kmers to use when counting the number "
                                    "of kmers in a file."
                           )
    arg_parser.add_argument( '-d', '--since', type = str,
                             help = "Beginning date to gather jobstats data from, "
                                    "if this argument is not provided, the current day will "
                                    "be used instead."
                           )

    # parse arguments
    args = arg_parser.parse_args()

    # Use today as the default date
    if not args.since:
        args.since = str( datetime.date.today() )

    # verify supplied arguments
    valid_args = verify_arguments( args )

    # if incorrect arguments were supplied, inform user and exit
    if valid_args != ErrorCodes.NO_ERROR:
        report_error( valid_args )
        sys.exit( 1 )

    # initialize variables 
    runner_obj   = SubprocessRunner()
    stats_parser = JobstatsParser()
    data_parser  = DataParser()
    kmer_parser  = KmerFileParser()
    command      = [ "jobstats", "-p", "-S %s" % args.since ]

    # call the command of the runner object
    command_result = runner_obj.invoke( command )

    # if the return is null, report that the return was non-zero, exit program
    if not command_result:
        report_error( ErrorCodes.NON_ZERO_EXIT_STATUS )
        sys.exit( 1 )

    # otherwise, set the data of the parser
    stats_parser.set_data( command_result )
    data_parser.set_data( args.job_data )
    kmer_parser.set_data( args.fasta_dir )

    # parse the data
    data_info     = DataInfo( data_parser.parse() )
    kmer_info     = KmerInfo( kmer_parser.parse(),
                              args.kmer_size, dir_name = args.fasta_dir
                            )
    jobstats_info = JobstatsInfo( stats_parser.parse() )

    jobstats_info.replace_job_name_with_filename( kmer_info.get_data() )

    # write to output file
    kmer_info.write( kmer_output )


class Parser( ABC ):
    def __init__( self ):
        pass

    @abstractmethod
    def parse( self ):
        pass

    @abstractmethod
    def set_data( self, new_data ):
        pass

class FileParser( Parser ):
    def __init__( self, filename = None ):
        self._filename = filename 

    def set_data( self, new_data ):
        self._filename = new_data 

class DataParser( FileParser ):
    def __init__( self, filename = None ):
        super().__init__( filename )

    def parse( self ):
        pass

class KmerFileParser( FileParser ):
    def __init__( self, filename = None ):
        super().__init__( filename )

    def parse( self ):
        with open( self._filename, 'r' ) as open_file:
            return self._parse_data( open_file )

    def _parse_data( self, file_obj ):
        out_dict = {}
        for line in file_obj:
            split_line = line.strip().split( '|' )
            cluster_name = split_line[ 0 ]
            job_id       = split_line[ 1 ]

            out_dict[ job_id ] = cluster_name 
        return out_dict
            
                           
class InfoClass:
    def __init__( self, data ):
        self._data = data

    def get_data( self ):
        return self._data 

    def write_to( self, filename ):
        with open( filename, 'w' ) as open_file:
            for item in self._data:
                open_file.write( str( item ) + '\n' )
        

class DataInfo( InfoClass ):
    def __init__( self, data ):
        super().__init__( data )
   
class JobstatsInfo( InfoClass ):
    class JobstatsData( Enum ):

        JOB_ID         = 0
        JOB_NAME       = 1
        JOB_REQ_MEM    = 2
        JOB_USED_MEM   = 3
        JOB_REQ_CPU    = 4
        JOB_USED_CPU   = 5
        JOB_TIME_LIMIT = 6
        JOB_ELAPSED    = 7
        JOB_STATE      = 8

    def __init__( self, data ):
        super().__init__( data )

    def replace_job_name_with_filename( self, kmer_dict ):
        stats_data = JobstatsInfo.JobstatsData
        for item in self._data:
            new_name = kmer_dict[ item[ stats_data.JOB_ID.value ] ]
            item[ stats_data.JOB_NAME.value ] = new_name
        print( self._data )

class KmerInfo( InfoClass ):
    def __init__( self, data, kmer_size, dir_name = None ):
        super().__init__( data )
        self._dir_name = dir_name
        self._kmer_size = kmer_size

        self._size_data = self._count_kmers()
        

    def _count_kmers( self ):
        out_dict = {}
        filenames = [ item for item in \
                      os.listdir( self._dir_name ) if item.endswith( '.fasta ' )
                    ]
        for current_file in filenames:
            out_dict[ current_file ] = self._count_kmers_in_file( current_file )
        return out_dict

    def _count_kmers_in_file( self, filename ):
        oligo_set = set()
        names, sequences = oligo.read_fasta_lists( filename )

        for item in sequences:
            oligo_set |= oligo.subset_lists_iter( item, self._kmer_size, 1 )
        return len( oligo_set )

    def write( self, outfile_name ):
        with open( outfile_name, 'w' ) as open_file:
            for key, value in self._size_data.items():
                open_file.write( "%s|%s\n" % ( key, value ) )
        
class JobstatsParser( Parser ):
    def __init__( self, data = None ):
        self._data = data

    def parse( self ):
        return_val = None

        if self._data:
            return_val = self._parse_data( self._data )
        return return_val

    def set_data( self, new_data ):
        self._data = new_data

    def _parse_data( self, str_data ):
        split_data = str_data.split( '\n' )
        split_data = [ item.strip().split( '|' ) for item in split_data if len( item ) > 0 ]
        return split_data
    
class Runner( ABC ):
    def __init__( self, runner = None ):
        self._runner       = runner
        self._error        = False
        self._last_command = ""

    @abstractmethod
    def invoke( self, command, args ):
        pass

    @abstractmethod
    def get_last_command( self ):
        pass

class SubprocessRunner( Runner ):
    def __init__( self ):
        super().__init__( runner = subprocess.check_output )

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

class ErrorCodes( Enum ):
    NO_ERROR                    = 0
    JOB_DATA_FILE_NOT_SUPPLIED  = 1
    JOB_DATA_FORMATTING_ERROR   = 2
    FASTA_DIR_NOT_FOUND         = 3
    INVALID_DATE_FORMAT         = 4
    NON_ZERO_EXIT_STATUS        = 5
    
def verify_arguments( args ):
    DATE_FORMAT = '%Y-%m-%d'
    if not args.job_data:
        return ErrorCodes.JOB_DATA_FILE_NOT_SUPPLIED
    try:
        datetime.datetime.strptime( args.since, DATE_FORMAT )
    except ValueError:
        return ErrorCodes.INVALID_DATE_FORMAT
    return ErrorCodes.NO_ERROR

def report_error( error_code ):
    errors = [ "No Error Found",
               "Job Data file not supplied",
               "Job Data file is formatted incorrectly",
               "Directory containing Fasta files either not found "
               "or is invalid",
               "Date supplied is not formatted correctly",
               "Jobstats returned non-zero exit status"
             ]

    print( "ERROR: %s. Program will exit..." % errors[ error_code.value ] )


if __name__ == '__main__':
    main()
