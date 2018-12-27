#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import subprocess                     # for getting jobstats data
from abc import ABC, abstractmethod   # for abstract base class
import protein_oligo_library as oligo # for getting kmer-counts


def main():
    pass
    # create parser object

    # add program options and arguments

    # parse arguments

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
