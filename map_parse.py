#!/usr/bin/env python3
import argparse

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

    # try to parse the map file
    try:
       # function: parse_map
        parse_map( args.map )
    # except error:
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

def parse_map( str: file_name ):
    """
        Opens, reads, and parses a file containing an epitope map.
        
        On error, raises the appropriate error.
     
        On successful operation, returns a dictionary containing
        epitope: list of items mapping
    """
    return 0

class InputFormatFileError( Exception ):
    pass
    

if __name__ == '__main__':
    main()
