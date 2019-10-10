#!/usr/bin/env python3
import sys
import argparse
import protein_oligo_library as oligo
import random
import os


def main():
    argp = argparse.ArgumentParser( description = 'Create some number of sub-fastas with '
                                                  'specified size '
                                                  'from a single fasta.'
                                  )
    argp.add_argument( '-i', '--input', type = str, help = "Name of the input fasta to sample. " )
    argp.add_argument( '-o', '--output', type = str, help = "Name of directory to write output files "
                       "to. If this directory already exists, the program will exit with an error."
                     )
    argp.add_argument( '--override_existing_dir', help = "Include if you want the script to continue "
                       "even if the specified output directory already exists.",
                       action = 'store_true'
                     )
    argp.add_argument( '-p', '--prefix', type = str, help = "The prefix for output"
                       "filenames. Each file that is written to the directory will "
                       "be written in the form prefix_XXXX, where XXXX is the 'ID' "
                       "of the sample file. This ID will start at 0 and go to the number "
                       "of output files produced. The number of digits necessary will be "
                       "determined.", default = "sample"
                     )
    argp.add_argument( '-n', '--num_sequences', help = "Then number of sequences "
                       "to include in each file.", type = str 
                     )
    argp.add_argument( '-f', '--num_samples', help = "The number of samples to draw from "
                       "the starting set of sequences."
                     )
    args = argp.parse_args()


    try:
        os.mkdir( args.output )
    except FileExistsError:
        if args.override_existing_dir:

            print( f"WARNING: Directory {args.output} exists, any files in the directory "
                   "whose names collide with the output files written by this script "
                   "will be overwritten."
                 )
        else:
            print( f"ERROR: Directory {args.output} exists, please remove the "
                   "existing directory or select another. "
                 )
            sys.exit( 1 )


if __name__ == '__main__':
    main()
