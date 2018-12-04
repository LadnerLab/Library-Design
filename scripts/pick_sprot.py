#!/usr/bin/env python3

import argparse                        # For parsing command-line arguments
import protein_oligo_library as oligo  # for operations on Fasta files


def main():
    arg_parser = argparse.ArgumentParser( description = "Parse representative output map to produce a sprot/trembl file" )

    arg_parser.add_argument( '-s', '--sprot', help = "Input sprot fasta to parse" )
    arg_parser.add_argument( '-t', '--trembl', help = "Input trembl file to parse" )
    arg_parser.add_argument( '-m', '--map_file', help = "Input map file to parse." )

    args = arg_parser.parse_args()

    out_sprot_name  = args.map_file + "_sprot"
    out_trembl_name = args.map_file + "_trembl"

    sprot_names, sprot_seqs   = oligo.read_fasta_lists( args.sprot )
    trembl_names, trembl_seqs = oligo.read_fasta_lists( args.trembl )

    in_sprot_seqs   = {}
    in_trembl_seqs  = {}
    out_sprot_seqs  = {}
    out_trembl_seqs = {}


    for index in range( len( sprot_names ) ):
        current_name = sprot_names[ index ]
        current_seq  = sprot_seqs[ index ]

        in_sprot_seqs[ current_name ] = current_seq

    for index in range( len( trembl_names ) ):
        current_name = trembl_names[ index ]
        current_seq  = trembl_seqs[ index ]

        in_trembl_seqs[ current_name ] = current_seq

    map_items = parse_map( args.map_file )
    print( map_items )


def parse_map( in_file ):
    out_list    = list()
    current_set = set()

    with open( in_file, 'r' ) as open_file:
        for line in open_file:
            split_line = line.strip().split( '\t' )
            current_set = set()
            current_set.update( split_line )
            out_list.append( current_set )
    return out_list 
            


if __name__ == '__main__':
    main()
