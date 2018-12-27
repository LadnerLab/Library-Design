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
        report_argument_error( valid_args )
        sys.exit( 1 )

    # initialize variables 

    # call the command of the runner object

    # if the return is null, report that the return was non-zero, exit program


    # otherwise, set the data of the parser

    # parse the data

    # write to output file

    
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
        super.__init__( subprocess.check_output )

    def invoke( self, command_arg_list ):
        self._last_command = command_arg_list

        try:
            value = self._runner( command_arg_list )
            # If we reach this point, runner returned exit status of zero
            self._error = False
            return value

        except subprocess.CalledProcessError:
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
    
def verify_arguments( args ):
    DATE_FORMAT = '%Y-%m-%d'
    if not args.job_data:
        return ErrorCodes.JOB_DATA_FILE_NOT_SUPPLIED
    try:
        datetime.datetime.strptime( args.since, DATE_FORMAT )
    except ValueError:
        return ErrorCodes.INVALID_DATE_FORMAT
    return ErrorCodes.NO_ERROR

def report_argument_error( error_code ):
    errors = [ "No Error Found",
               "Job Data file not supplied",
               "Job Data file is formatted incorrectly",
               "Directory containing Fasta files either not found "
               "or is invalid",
               "Date supplied is not formatted correctly"
             ]

    print( "ERROR: %s. Program will exit..." % errors[ error_code.value ] )


if __name__ == '__main__':
    main()
