#!/usr/bin/env python3
import sys

def main():

    nodes_dmp   = sys.argv[ 1 ]
    lineage_dmp = sys.argv[ 2 ]

    nodes_dict = parse_nodes( nodes_dmp )


def parse_nodes( file_name ):
    READ_ONLY_FLAG = "r"
    nodes_file     =  open( file_name,
                            READ_ONLY_FLAG
                          )
    nodes_dict = {}

    for line in nodes_file:
        continue

    nodes_file.close()

    return nodes_dict

if __name__ == '__main__':
    main()
