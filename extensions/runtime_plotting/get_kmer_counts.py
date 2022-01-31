#!/usr/bin/env python3
import protein_oligo_library as oligo
from glob import glob
import os
import sys
import argparse

def main():
    argp = argparse.ArgumentParser( "Count the number of kmers in each file in a directory " 
                                    "for each directory in a starting directory."
                                  )

    argp.add_argument( "-i", "--in_dir", type = str,
                       help = "The input directory to scan for "
                              "directories containing fasta files."
                     )

    argp.add_argument( '-k', '--kmer_size',
                       type = int,
                       help = "The kmer size to consider "
                              "when counting kmers."
                     )

    argp.add_argument( '-o', '--output', type = str,
                       help = "The name of the file to write "
                              "output to. This will be a tab-delimited "
                              "file with no header. The first entry in a "
                              "tab-delimited column will be the path of the file "
                              "relative to the directory in which this script was invoked. "
                              "The second will be the number of kmers in the file."
                     )

    argp.add_argument( '-t', '--total', action = 'store_true',
                       help = "Include this flag if the output should be the "
                              "total number of kmers in each file, not the "
                              "number of unique kmers."
                     )

    args = argp.parse_args()
    in_dir = args.in_dir
    k = args.kmer_size
    out_file = args.output
    output_total = args.total
    

    in_files = glob( f'{in_dir}/*')

    with open( out_file, 'w' ) as out_f:
        count_kmers_in_dir( out_f,
                            in_files,
                            k,
                            count_total = output_total
                          )

def count_kmers_in_dir( out_fh, list_of_files, k, count_total = False ):
    
    for fname in list_of_files:
        if os.path.isdir( fname ):
            count_kmers_in_dir( out_fh,
                                glob( f'{fname}/*' ), k,
                                count_total = count_total
                              )
        else:
            kmer_count = count_kmers_in_file( fname, k,
                                              count_total = count_total
                                            )

            out_fh.write( f'{fname}\t{kmer_count}\n' )

def count_kmers_in_file( filename, k, count_total = False ):
    names, sequences = oligo.read_fasta_lists( filename )
    kmers = set()
    ret_val = 0

    for seq in sequences:
        kmer_set = oligo.subset_lists_iter( seq, k, 1 )

        if count_total:
            ret_val += len( kmer_set )
        else:
            kmers |= kmer_set
            ret_val = len( kmers )

    return ret_val

if __name__ == '__main__':
    main()
