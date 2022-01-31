#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo

def main():
    original = sys.argv[ 1 ]
    clusters = sys.argv[ 2 ]

    orig_dict = {}

    orig_names, orig_seqs = oligo.read_fasta_lists( original )
    clusters_names, clusters_seqs = oligo.read_fasta_lists( clusters )

    assert len( orig_names ) == len( orig_seqs )
    assert len( clusters_names ) == len( clusters_seqs )
    for current in range( len( orig_names ) ):
        orig_dict[ orig_names[ current ] ] = orig_seqs[ current ]

    successful = 0
    for current in range( len( clusters_names ) ):
        name = clusters_names[ current ]
        seq  = clusters_seqs[ current ]

        if name in orig_dict:
            assert orig_dict[ name ] == seq
            successful += 1






if __name__ == '__main__':
    main()
