#!/usr/bin/env python3
import argparse
import sys

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Script to parse information given by validate_design's epitope map output" )

    arg_parser.add_argument( '-m', '--map', help = "Name of file containing epitope map" )
    arg_parser.add_argument( '-o', '--output', default = "parsed_map",
                             help = (
                                      "Identifier to preface output files."
                                      " Output will be written in a tab-delimited format"
                                    )
                           )
    arg_parser.add_argument( '-t', '--tax_db',
                             help = "Name of file containing mappings of taxonomic id -> rank data"
                           )
    arg_parser.add_argument( '-v', '--verbose',
                             help = "Flag to add if output should be written to STDOUT"
                           )
    arg_parser.add_argument( '-g', '--gap_file',
                             help = ( "File containing mappings of taxid->rank for "
                                      "use in filling gaps."
                                    )
                           )
    arg_parser.add_argument( '-r', '--reference',
                             help = "Fasta reference dataset used to create library"
                           )

    args = arg_parser.parse_args()

    oligo_centric_table    = None
    sequence_centric_table = None
    species_centric_table  = None

    map_dict = None
    gap_dict = parse_gaps( args.gap_file )

