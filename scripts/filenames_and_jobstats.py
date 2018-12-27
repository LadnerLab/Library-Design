#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import subprocess                     # for getting jobstats data
from abc import ABC, abstractmethod   # for abstract base class
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

    # verify supplied arguments

    # if incorrect arguments were supplied, inform user and exit

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

if __name__ == '__main__':
    main()
