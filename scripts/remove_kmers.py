#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import sys                            # for exiting upon failure
from enum import Enum                 # for error reporting
import protein_oligo_library as oligo # for operations on oligos/fasta files


def main():

    parser = argparse.ArgumentParser( description = "Given design and oligo-centric "
                                                    "table, remove kmers from a design "
                                                    "that are found in multiple "
                                                    "of certain categories."
                                    )
    parser.add_argument( '-l', '--lower_bound', type = int, default = 0,
                         help = "Lower-bound (inclusive) for a value to be "
                                "considered for output. If a value is less "
                                "than this, an item is out for consideration. "
                                "Note that an item must be greater than lower "
                                "bound and less than upper bound to be considered."
                                "Also note that lower_bound <= upper_bound."
                                "To include items that are exactly a value, set "
                                "upper and lower bound to the same value."

                     )

    parser.add_argument( '-u', '--upper_bound', type = int, default = 1,
                         help = "Upper-bound (inclusive) for a value to be "
                                "considered for output. If a value is greather  "
                                "than this, it will not be considered for output. "
                                "Note that lower_bound <= upper_bound."
                                "To include items that are exactly a value, set "
                                "upper and lower bound to the same value."
                     )

    parser.add_argument( '-t', '--table', type = str,
                         help = "Oligo-centric table ( as created by map_parse ) "
                                "to parse. Values in the table, "
                                "as specified by command-line arguments, "
                                "are parsed and the oligos that 'fit' the supplied "
                                "parameters are chosen."
                              

                     )

    parser.add_argument( '-d', '--design', type = str,
                         help = "Name of the fasta file containing "
                                "the designed oligos."
                     )

    # parse options and arguments
    args = parser.parse_args()

    UPPER_BOUND = args.upper_bound
    LOWER_BOUND = args.lower_bound

    validate_args( args )

    # read the oligos from the library design into a dict
        # function: oligo.sequence_dict_from_file

    # try to read and parse the input file
        # function: read_oligo_table

    # on exception, report error to the user

    # assume valid input data

    # 

class CommandArgError( Enum ):
    NO_ERROR                  = 0,
    TABLE_NOT_SUPPLIED_ERROR  = 1,
    DESIGN_NOT_SUPPLIED_ERROR = 2,
    INVALID_BOUNDS_ERROR      = 3,
    TABLE_FILE_IO_ERROR       = 4,
    DESIGN_FILE_IO_ERROR      = 5
    
def validate_args( args ):
    """
        Ensures that args supplied by user are valid.
    
        :pre: upper_bound < lower_bound
        :pre: table was supplied by the command line
        :pre: design was supplied by the command line
    
        :note: Does not check that the files table and design actually exist,
               only that they were provided by the user.
    
        :post: Returns CommandArgError.NO_ERROR on success, 
               on failure reports the appropriate error. 
               
    """

    # check that upper_bounds is greater than lower_bounds
    if args.upper_bound < args.lower_bound:
        return CommandArgError.INVALID_BOUNDS_ERROR

    # check that table was supplied
    elif args.table is None:
        return CommandArgError.TABLE_NOT_SUPPLIED_ERROR

    # check that design was supplied
    elif args.design is None:
        return CommandArgError.DESIGN_NOT_SUPPLIED_ERROR
    
    return CommandArgError.NO_ERROR

if __name__ == '__main__':
    main()
