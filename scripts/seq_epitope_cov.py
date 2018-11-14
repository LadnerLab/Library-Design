#!/usr/bin/env python3
import argparse                        # for parsing command-line arguments
import sys                             

import protein_oligo_library as oligo # for operations on sequences

def main():
    parser = argparse.ArgumentParser( description = "Calculate sequence-level epitope coverage for "
                                                    "some library design."
                                    )

    parser.add_argument( '-d', '--design', help = "Name of design file to compare against "
                                                  "reference file."
                       )
    parser.add_argument( '-r', '--reference', help = "Reference dataset from which 'design' was "
                                                     "created, script calculates epitope coverage "
                                                     "for each sequence in the reference."
                       )
    parser.add_argument( '-o', '--output', help = "Name of file to write output to. " )
    parser.add_argument( '--delimiter', default = '|', help = "Delimiter to place between sequence name "
                                                              "and its percent epitope coverage."
                       )
    parser.add_argument( '-e', '--epitope_size', type = int,
                         help = "Integer size of epitopes to use for coverage "
                                "information gathering."
                       )

    args = parser.parse_args()

    ref_names, ref_sequences  = oligo.read_fasta_lists( args.reference )
    design_names, design_seqs = oligo.read_fasta_lists( args.design )

    header = "Sequence Name%sPercent Epitopes (%d-mers) Covered in design\n" % (
             args.delimiter, args.epitope_size
             )

    out_strings  = list()
    design_kmers = set()
    ref_dict     = create_kmer_dict( ref_names, ref_sequences, args.epitope_size )

    for item in design_seqs:
        design_kmers |= oligo.subset_lists_iter( item, args.epitope_size, 1 )

    for current_item in ref_dict:
        current_kmers = ref_dict[ current_item ]
        perc_cov = len( current_kmers & design_kmers ) / len( current_kmers )

        out_strings.append( "%s%s%f\n" % ( current_item, args.delimiter, perc_cov ) )

    with open( args.output, 'w' ) as out_file:
        out_file.write( "%s\n" % header )
        for item in out_strings:
            out_file.write( item )
        

def create_kmer_dict( names, sequences, epitope_size ):
    out_dict = {}

    for index in len( range( names ) ):
        current_name     = names[ index ]
        current_sequence = sequences[ index ]

        kmers = oligo.subset_lists_iter( current_sequence,
                                         epitope_size,
                                         1
                                       )
        out_dict[ current_name ] = kmers 

    return out_dict

        
    
                          


if __name__ == '__main__':
    main()
