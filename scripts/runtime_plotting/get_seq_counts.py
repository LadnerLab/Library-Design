#!/usr/bin/env python3
import protein_oligo_library as oligo
from glob import glob
import os
import sys

def main():

    if len( sys.argv ) != 4:
        print( "USAGE: get_kmer_counts in_dir k out_file")
        sys.exit( 1 )

    in_dir = sys.argv[ 1 ]
    k = int( sys.argv[ 2 ] )
    out_file = sys.argv[ 3 ]

    in_files = glob( f'{in_dir}/*')

    with open( out_file, 'w' ) as out_f:
        count_kmers_in_dir( out_f, in_files, k )

def count_kmers_in_dir( out_fh, list_of_files, k ):
    
    for fname in list_of_files:
        if os.path.isdir( fname ):
            count_kmers_in_dir( out_fh, glob( f'{fname}/*' ), k )
        else:
            kmer_count = count_kmers_in_file( fname, k )

            out_fh.write( f'{fname}\t{kmer_count}\n' )

def count_kmers_in_file( filename, k ):
    names, sequences = oligo.read_fasta_lists( filename )
    kmers = set()

    for seq in sequences:
        kmers |= oligo.subset_lists_iter( seq, k, 1 )
    return len( kmers )

if __name__ == '__main__':
    main()
