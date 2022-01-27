#!/usr/bin/env python3

import argparse
import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Verify a map that has been created by the one_hundred_reps script, checks "
                                                        "to ensure that all collapsed sequences were collapsed correctly"
                                        )
    arg_parser.add_argument( '-u', '--uncollapsed_input', help = "FASTA containing uncollapsed sequences" )
    arg_parser.add_argument( '-m', '--map_file', help = "Map file produced by one_hundred_reps containing "
                                                        "sequence to collapsed sequences mapping"
                           )

    args = arg_parser.parse_args()

    original_seq_dict    = fasta_to_dict( args.uncollapsed_input )
    collapsed_names_dict = parse_map( args.map_file )

    original_len   = len( original_seq_dict.keys() )
    collapsed_len  = len( collapsed_names_dict.keys() )

    print( "Original number of sequences:  %d." %  original_len )
    print( "Collapsed number of sequences: %d." %  collapsed_len )

    print( "Number of sequences who were collapsed (expected): %d." % ( original_len - collapsed_len ) )
    print( "Number of sequences who were collapsed (actual):   %d." % sum( [ len( item ) for item in collapsed_names_dict.values() ] ) )
    
           

    for collapsed_under, collapsed in collapsed_names_dict.items():
        for current in collapsed:
            if original_seq_dict[ current ] not in original_seq_dict[ collapsed_under ]:
                print( original_seq_dict[ current ], original_seq_dict[ collapsed_under ] )



def fasta_to_dict( filename ):
    names, sequences = oligo.read_fasta_lists( filename )
    out_dict = {}

    for index, name in enumerate( names ):
        out_dict[ name ] = sequences[ index ]
    return out_dict

def parse_map( map_filename ):
    out_dict = {}
    with open( map_filename, 'r' ) as open_file:
        for line in open_file:
            names = line.strip().split( '\t' )
            new_key = names[ 0 ].strip().replace( '>', '' )

            out_dict[ new_key ] = list()
            for name in names[1::]:
                out_dict[ new_key ].append( name.replace( '>', '' ) )
    return out_dict
                

if __name__ == '__main__':
    main()
