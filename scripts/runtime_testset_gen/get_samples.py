#!/usr/bin/env python3
import argparse
import os

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
    argp.add_argument( '-q', '--quartiles', help = "Comma-separated quartiles to grab. Quintile 100 "
                       "will grab the largest, and quintile 0 the smallest."
                     )
    argp.add_argument( '--override_existing_dir', help = "Include if you want the script to continue "
                       "even if the specified output directory already exists.",
                       action = 'store_true'
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
