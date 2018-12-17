#!/usr/bin/env python3


import protein_oligo_library as oligo
import sys


def main():
    in_file = sys.argv[ 1 ]

    names, sequences = oligo.read_fasta_lists( in_file )

    print( "Number of oligos in original design: %d" % len( names ) )
    print( "Number of unique oligos:             %d" % len( set( sequences ) ) )

    names, sequences = oligo.get_unique_sequences( names, sequences )

    oligo.write_fastas( names, sequences, in_file + "_unique" )




if __name__ == '__main__':
    main()


