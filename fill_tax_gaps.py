#!/usr/bin/env python3
import sys

def main():

    nodes_dmp   = sys.argv[ 1 ]
    lineage_dmp = sys.argv[ 2 ]
    out_file    = sys.argv[ 3 ]

    nodes_dict = parse_nodes( nodes_dmp )


def parse_nodes( file_name ):
    READ_ONLY_FLAG = 'r'
    DELIMITER_CHAR = '|'
    nodes_file     =  open( file_name,
                            READ_ONLY_FLAG
                          )
    nodes_dict = {}

    current_id        = 0
    current_parent_id = 0
    current_rank      = ""

    for line in nodes_file:
        split_line = [ item.strip() \
                       for item in \
                       line.split( DELIMITER_CHAR )
                     ]  

        current_id        = int( split_line[ 0 ] )
        current_parent_id = int( split_line[ 1 ] )
        current_rank      = split_line[ 2 ]



    nodes_file.close()

    return nodes_dict

if __name__ == '__main__':
    main()
