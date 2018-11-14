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

    for item in design_seqs:
        design_kmers |= oligo.subset_lists_iter( item, args.epitope_size, 1 )

    for index in range( len( ref_names ) ):
        current_name  = ref_names[ index ]
        current_seq   = ref_sequences[ index ]
        current_kmers = oligo.subset_lists_iter( current_seq, args.epitope_size, 1 )

        if len( current_kmers ) > 0:
            perc_cov      = len( current_kmers & design_kmers ) / len( current_kmers )
        else:
            perc_cov = 4

        out_strings.append( "%s%s%f\n" % ( current_name, args.delimiter, perc_cov ) )

    with open( args.output, 'w' ) as out_file:
        out_file.write( "%s\n" % header )
        for item in out_strings:
            out_file.write( item )

if __name__ == '__main__':
    main()
