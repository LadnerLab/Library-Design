#!/usr/bin/env python3
import argparse
import sys
import protein_oligo_library as oligo
import re


def main():
    argparser = argparse.ArgumentParser( description = "Retrieve all sequences in a file that "
                                                       "have contain a specified taxonomic id."
                                       )
    argparser.add_argument( '-i', '--ids', help = "Name of file containing one "
                                                  "id per line"
                          )
    argparser.add_argument( '-f', '--fasta', help = "Name of fasta file whose sequences "
                                                    "will be searched for the specified ids."
                          )
    args = argparser.parse_args()

    seq_dict        = oligo.sequence_dict_from_file( args.fasta )
    ids_of_interest = get_id_set( args.ids )
    taxid_dict      = get_tax_ids( seq_dict )

    query_pattern = r'OXX=[0-9]* *,[0-9]* *,[0-9]* *,[0-9]* *'
            tax_ids = re.search( query_pattern, query ).group()

def get_id_set( filename ):
    out_set = set()
    with open( filename, 'r' ) as open_file:
        for line in open_file:
            out_set.add( line.strip() )
    return out_set
            

if __name__ == '__main__':
    main()
