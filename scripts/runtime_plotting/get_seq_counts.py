#!/usr/bin/env python3
import protein_oligo_library as oligo
from glob import glob
import os
import sys

def main():

    if len( sys.argv ) != 3:
        print( f"USAGE: {sys.argv[0]} in_dir out_file")
        sys.exit( 1 )

    in_dir = sys.argv[ 1 ]
    out_file = sys.argv[ 2 ]

    in_files = glob( f'{in_dir}/*')

    with open( out_file, 'w' ) as out_f:
        count_sequences_in_dir( out_f, in_files )

def count_sequences_in_dir( out_fh, list_of_files ):
    
    for fname in list_of_files:
        if os.path.isdir( fname ):
            count_sequences_in_dir( out_fh, glob( f'{fname}/*' ) )
        else:
            sequence_count = count_sequences_in_file( fname )

            out_fh.write( f'{fname}\t{sequence_count}\n' )

def count_sequences_in_file( filename ):
    names, sequences = oligo.read_fasta_lists( filename )
    return len( sequences )

if __name__ == '__main__':
    main()
