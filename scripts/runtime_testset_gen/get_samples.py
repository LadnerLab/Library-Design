#!/usr/bin/env python3
import argparse
import os
import protein_oligo_library as oligo
import numpy as np
import random
from shutil import copyfile

def main():
    argp = argparse.ArgumentParser( description = "" )

    argp.add_argument( '-i', '--in_dir', help = "Name of the input directory to grab "
                       "files from."
                     )
    argp.add_argument( '-o', '--output_dir', help = "Name of the input directory to write "
                       "output files to."
                     )
    argp.add_argument( '-k', '--kmer_size', help = "The kmer size to use.", type = int,
                       default = 7
                     )
    argp.add_argument( '-q', '--percentiles', help = "Comma-separated percentiles to grab. Quantiles 100 "
                       "will grab the largest, and quintile 0 the smallest."
                     )
    argp.add_argument( '--override_existing_dir', help = "Include if you want the script to continue "
                       "even if the specified output directory already exists.",
                       action = 'store_true'
                     )
    argp.add_argument( '-t', '--total_files', help = "The total number of files to generate. ",
                       type = int
                       )

    args = argp.parse_args()
    quantiles = sorted( list( map( lambda x: float( x ) / 100.0,
                                   args.percentiles.split( ',' )
                                 )
                            )
                      )

    try:
        os.mkdir( args.output_dir )
    except FileExistsError:
        if args.override_existing_dir:

            print( f"WARNING: Directory {args.output_dir} exists, any files in the directory "
                   "whose names collide with the output files written by this script "
                   "will be overwritten."
                 )
        else:
            print( f"ERROR: Directory {args.output_dir} exists, please remove the "
                   "existing directory or select another. "
                 )
            sys.exit( 1 )

    fname_w_counts = get_kmer_sizes( args.in_dir, args.kmer_size )

    files_sorted = sorted( fname_w_counts, key = lambda x: x[ 1 ] )
    files_only_counts = [ x[ 1 ] for x in files_sorted ]
    file_quantiles = np.quantile( files_only_counts,
                                   quantiles,
                                   interpolation = 'nearest'
                                 ) 
    outputs = list()
    for index, value in enumerate( file_quantiles[:-1] ):
        start_idx = list( files_only_counts ).index( value )
        end_idx   = list( files_only_counts ).index( file_quantiles[ index + 1 ] )
        selected = random.sample( files_sorted[ start_idx : end_idx ],
                                  args.total_files // ( len( quantiles ) - 1)
                                )
        outputs.append( selected )

    # always grab the largest and smallest one
    outputs.append( [ files_sorted[ 0 ], files_sorted[ -1 ] ] )

    for out_group in outputs:
        for f_out in out_group:
            out_f = f_out[ 0 ]
            copyfile( f'{args.in_dir}/{out_f}', f'{args.output_dir}/{out_f}' )


def get_kmer_sizes( input_dir, k ):
    output = list() 
    for f_in in os.listdir( input_dir ):
        fname = f'{input_dir}/{f_in}'

        names, seqs = oligo.read_fasta_lists( fname )
        output.append( ( f_in, seqs_to_kmers( seqs, k ) ) )

    return output

def seqs_to_kmers( seqs, k ):
    output = set()
    for seq in seqs:
        output |= oligo.subset_lists_iter( seq, k, 1 )
    return len( output )


if __name__ == '__main__':
    main()
