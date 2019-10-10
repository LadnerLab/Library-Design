#!/usr/bin/env python3
import sys
import argparse
import protein_oligo_library as oligo
import random
import math
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
    argp.add_argument( '-s', '--sample_size', help = "Then number of sequences "
                       "to include in each file.", type = int 
                     )
    argp.add_argument( '-n', '--num_samples', help = "The number of samples to draw from "
                       "the starting set of sequences.", type = int
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

    input_seqs = set( parse_fasta( args.input ) )

    samples = sample_seqs( input_seqs,
                           args.num_samples,
                           args.sample_size
                         )

    write_samples( args.output, args.prefix,
                   samples
                 )

def write_samples( output_dir, file_prefix,
                   samples
                 ):
    num_digits = math.floor( math.log10( len( samples ) ) ) + 1

    digit_str = lambda x: str( x ).zfill( num_digits )

    for index, sample in enumerate( samples ):
        filename = f'{output_dir}/{file_prefix}_{digit_str( index )}'

        write_sample( filename, sample )


def sample_seqs( sequences,
                 num_samples,
                 samplesize,
                 replacement = False
               ):

    return_sets = list()
    # we don't want to sample more than we have
    total_sample_size = min( len( sequences ),
                             num_samples * samplesize
                           )

    if replacement:
        # TODO
        pass
    else:
        total_sample = random.sample( sequences, total_sample_size )

    for index in range( 0, total_sample_size, samplesize ):
        sample = set( total_sample[ index : index + samplesize ] )
        return_sets.append( sample )

    # fewer sequences than samples? 
    if total_sample_size == len( sequences ):
        remainder = len( total_sample ) % samplesize
        # grab the last sequences
        last_sample = set( total_sample[ -1 : -( remainder + 1 ) : -1 ] )
        return_sets.append( last_sample )

    return return_sets

    

def parse_fasta( fname ):
    names, sequences = oligo.read_fasta_lists( fname ) 

    return [ Sequence( a, b ) for a, b in zip( names, sequences ) ]

class Sequence:
    def __init__( self, name, seq ):
        self._name = name
        self._seq = seq

    def __hash__( self ):
        return hash( self._seq )

    def __eq__( self, other ):
        return self.z_seq == other._seq

    def __ne__( self, other ):
        return not ( self == other )

    def get_name( self ):
        return self._name

    def get_seq( self ):
        return self._seq

    def set_seq( self, seq ):
        self._seq = seq

    def set_name( self, name ):
        self._name = name

if __name__ == '__main__':
    main()
