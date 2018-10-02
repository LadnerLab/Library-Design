#!/usr/bin/env python3
import argparse
import sys

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Script to parse information given by validate_design's epitope map output" )

    arg_parser.add_argument( '-m', '--map', help = "Name of file containing epitope map" )
    arg_parser.add_argument( '-o', '--output', default = "parsed_map.txt",
                             help = (
                                      "Name of file to write parsed output to."
                                      " Output will be written in a tab-delimited format"
                                    )
                           )
    arg_parser.add_argument( '-t', '--tax_db',
                             help = "Name of file containing mappings of taxonomic id -> rank data"
                           )
    arg_parser.add_argument( '-v', '--verbose',
                             help = "Flag to add if output should be written to STDOUT"
                           )

    args = arg_parser.parse_args()

    map_dict = None
    try:
        map_dict = parse_map( args.map )
    except ( IOError, OSError, TypeError ):
        print( "ERROR: An IO exception occurred when trying to "
               "open and parse map file"
             )
        sys.exit( 1 )
    except InputFormatFileError:
        print( "ERORR: The input file provided "
               "is not formatted correctly and cannot be "
               "parsed by this program."
             )
        sys.exit( 1 )

    #
    print( map_dict )

def parse_map( file_name ):
    """
        Opens, reads, and parses a file containing an epitope map.
        
        On error, raises the appropriate error.
     
        On successful operation, returns a dictionary containing
        epitope: list of items mapping
    """
    READ_FLAG  = 'r'
    TAB_CHAR   = '\t'
    COMMA_CHAR = ','

    file_dict = {}
    open_file = None

    split_line   = ""
    new_dict_key = ""
    new_dict_val = ""

    # try to open the input file
    try:
        open_file = open( file_name, READ_FLAG )
        # loop through each line in file
        for line in open_file:
            # split the line on the tab character
            try:
                split_line = line.split( TAB_CHAR )
                new_dict_key = remove_loc_markers( split_line[ 0 ] )
                new_dict_val = split_line[ 1 ]
            except ( IndexError, Exception ):
                raise InputFormatFileError

            if len( new_dict_key ) == 0 \
               or len( new_dict_val ) == 0:
                raise InputFormatFileError
            
            # String is formatted correctly 
            else:
                if new_dict_key in file_dict.keys():
                    file_dict[ new_dict_key ] = file_dict[ new_dict_key ] + \
                                                new_dict_val.strip().split( COMMA_CHAR )
                else:
                    file_dict[ new_dict_key ] = new_dict_val.strip().split( COMMA_CHAR )
                    
    except:
        raise

    open_file.close()

    return file_dict

def remove_loc_markers( input_str ):
    """
        Removes the location markers from a name in the epitope map file

        Input is the name of a sequences followed by int_int
     
        e.g. Seq1_23_45
    """
    UNDERSCORE_CHAR = '_'

    split_str = input_str.split( UNDERSCORE_CHAR )

    if len( split_str ) < 3:
        raise InputFormatFileError

    split_str.pop()
    split_str.pop()

    return '_'.join( split_str )

class InputFormatFileError( Exception ):
    pass
    

if __name__ == '__main__':
    main()
