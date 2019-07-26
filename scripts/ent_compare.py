#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "" )

    arg_parser.add_argument( '-k', '--kmer_size', help = "The kmer size to use for evaluations." )
    arg_parser.add_argument( '-d', '--designed', help = "Fasta file containing designed sequences." )
    arg_parser.add_argument( '-c', '--check', help = "Fasta file containing sequence to check coverage." )
    arg_parser.add_arguemnt( '--create_table', action = 'store_true' )

    args = arg_parse.parse_args()

    # if a table is created choose to do that
    if args.create_table:
        create_table( args )


def create_table( arguments ):
    pass

if __name__ == '__main__':
    main()
