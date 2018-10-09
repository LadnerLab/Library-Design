#!/usr/bin/env python3
import sys

import protein_oligo_library as oligo

def main():

    nodes_dmp   = sys.argv[ 1 ]
    out_file    = sys.argv[ 2 ]

    nodes_dict = parse_dict( nodes_dmp )

    write_dict( nodes_dict, out_file )

def parse_dict( nodes_file ):
    READ_ONLY_FLAG = 'r'
    DELIMITER_CHAR = '|'
    nodes_file     =  open( nodes_file,
                            READ_ONLY_FLAG
                          )
    nodes_dict = {}

    current_id        = 0
    current_rank      = ""

    for line in nodes_file:
        
        split_line = [ item.strip() \
                       for item in \
                       line.split( DELIMITER_CHAR )
                     ]  

        current_id        = int( split_line[ 0 ] )
        current_rank      = split_line[ 2 ]

        nodes_dict[ current_id ] = current_rank.upper().replace( ' ', '_' )

    nodes_file.close()

    return nodes_dict

def write_dict( parsed_dict, out_file ):
    WRITE_FLAG = 'w'

    out_str = ""

    for current in parsed_dict.keys():
        out_str += "%s|%s\n" % ( current, parsed_dict[ current ] )
    
    open_file = open( out_file, WRITE_FLAG )

    open_file.write( out_str )
    open_file.close()

if __name__ == '__main__':
    main()
