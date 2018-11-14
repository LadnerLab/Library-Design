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
                          


if __name__ == '__main__':
    main()
