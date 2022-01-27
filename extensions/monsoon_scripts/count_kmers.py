#!/usr/bin/env python3

import sys
import protein_oligo_library as oligo

def main():
    names, sequences = oligo.read_fasta_lists( sys.argv[ 1 ] ) 

    kmer_set = set()
    for seq in sequences:
        kmer_set |= oligo.subset_lists_iter( seq, 9, 1 )
    print( "%s|%d" % ( sys.argv[ 1 ], len( kmer_set ) ) )


if __name__ == '__main__':
    main()
