#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse
import random


def main():
    argp = argparse.ArgumentParser( description = "Randomly generate a fasta "
                                                  "containing AA sequences."
                                  )
    argp.add_argument( '-f', '--filename' )
    argp.add_argument( '-l', '--length', type = int )
    argp.add_argument( '-g', '--generate',
                       type = int
                     )
    argp.add_argument( '-p', '--prefix', help = "The prefix for each "
                                                "generated sequence. ",
                       default = "sequence"
                     )

    args = argp.parse_args()

    names, sequences = list(), list()

    for index in range( args.generate ):
        name = f'{args.prefix}_{index}'
        sequence = generate_aa_sequence( args.length )

        names.append( name ) 
        sequences.append( sequence ) 

    oligo.write_fastas( names, sequences, output_name = args.filename )


def generate_aa_sequence( length ):
    candidate_aa = ['G', 'K', 'N', 'R', 'C', 'H',
                    'I', 'L', 'W', 'E', 'B', 'D',
                    'A', 'V', 'U', 'S', 'P', 'F',
                    'Q', 'Z', 'Y', 'T', 'M'
                   ]
    return ''.join( random.choices( candidate_aa, k = length ) )

if __name__ == '__main__':
    main()
