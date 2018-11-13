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
    parser.add_argument( '-c', '--consideration', action = 'append',
                         help = "Values to consider when removing oligos. "
                                "( adding oligos to the output ). For an "
                                "oligo to be included, it must have at least "
                                "lower_bound and at most upper_bound of this "
                                "and all items added with this flag. "
                                "Including this flag multiple times with different "
                                "values results in more items being checked."
                                "Valid options include: "
                                "num_seqs, num_species, num_genera, "
                                "num_families."
                       )

    # parse options and arguments
    args = parser.parse_args()

    UPPER_BOUND = args.upper_bound
    LOWER_BOUND = args.lower_bound

    args_result = validate_args( args )

    if args_result != CommandArgError.NO_ERROR:
        report_error( args_result )
        sys.exit( 1 )

    # read the oligos from the library design into a dict
    design_dict = parse_fasta_file( args.design )

    if not design_dict:
        report_error( CommandArgError.DESIGN_FILE_IO_ERROR )

        sys.exit( 1 )


    # try to read and parse the input file
    table_dict = parse_oligo_table( args.table )

    if not table_dict:
        report_error( CommandArgError.TABLE_FILE_IO_ERROR )

        sys.exit( 1 )

    # assume valid input data

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

def parse_fasta_file( file_name ):
    """
        Reads and parses a fasta file, 
        located at file_name. 
        Upon successful operation, returns 
        a dictionary containing name: sequence 
        mappings parsed from file_name. Upon IOError,
        returns None.
    """

    return_dict = None

    try:
        # read the files stored in file_name into a dict
        return_dict = oligo.sequence_dict_from_file( file_name )
    except IOError:
        return_dict = None
    finally:
        return return_dict

def report_error( error_code ):
    """
        Reports an error based upon the CommandArgError value supplied.
    """
    error_strings = [ "No error",
                      "An oligo table was not supplied",
                      "The design fasta file was not supplied",
                      "Invalid bounds: upper_bound must not be less than lower_bound",
                      "The supplied table file was not found, "
                      "or another error occurred when opening the file",
                      "The supplied design fasta file was either not found "
                      "or another error occurred when trying to access it"
                      ]

    print( "ERROR: %s. Please resolve this "
           "issue and restart the script. " % error_strings[ error_code.value[ 0 ] ]
         )

def parse_oligo_table( filename, species_covered = False ):
    """
        Parses the oligo-centric table file found at filename.
        
        :post: Upon success, returns a dictionary containing mappings of:
               oligo: [ num_seqs, num_species, num_genera, num_families ]
               If species_covered is set to True, then a list of the names of the 
               species covered by this oligo will be 4th element in each value.

        :post: Upon failure (Due to IO or improper format), 
               returns None
    """

    DELIMITER_CHAR = '\t'

    return_dict = {}

    try:
        open_file = open( filename, 'r' )

        # first line contains headers, skip it
        for current_entry in list( open_file )[ 1:: ]:
            split_line = current_entry.strip( '\n' ).split( DELIMITER_CHAR )

            oligo_name = split_line[ 0 ].strip()

            return_dict[ oligo_name ] = get_items_from_entry( split_line, species_covered )

    # upon error, set return_val to None
    except( IOError, ValueError ):
        return_dict = None

    # finally
    finally:
        return return_dict

def get_items_from_entry( line_list, taxons_covered ):
    """
        Gets the appropriate items from line_list, casts them as integer.
        If a list does not contain the same number of entries as all of the others,
        raises ValueError as the table is formatted incorrectly. 
 
        :pre: line_list contains one line from an oligo_table that is nto the header
        :pre: species_covered is either True or False.

        :post: returns a list, with entries being:
               [ num_sequences, num_species, num_genera, num_families ]
               If species_covered is True, the last item in the list is a list of all of the 
               species some entry covered
    """

    DELIMITER_CHAR = ','

    return_val = None

    # raises ValueError if this does not work, thus our data is validated here
    try:
        name, num_seqs, num_species, \
        num_genera, num_fam, \
        species_cov, genera_cov, family_cov = line_list

    except ValueError:
        return None

    return_val = [ int( num_seqs ), int( num_species ), int( num_genera ),
                   int( num_fam )
                 ]

    if taxons_covered:
        return_val.append( species_cov.split( ',' ) )
        return_val.append( genera_cov.split( ',' ) )
        return_val.append( family_cov.split( ',' ) )

    return retun_val
        

if __name__ == '__main__':
    main()
