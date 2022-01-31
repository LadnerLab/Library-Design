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
    target_names    = get_desired_names( taxid_dict, ids_of_interest )

    print( "Number of sequences containing ids: %d." % len( target_names ) )

    for name in target_names:
        print( ">%s\n%s" % ( name, seq_dict[ name ] ) )

def get_desired_names( id_dict, id_set ):
    out_names = set()

    for id in id_set:
        out_names |= id_dict[ id ]
    return out_names
    
def get_tax_ids( seq_dict ):
    id_dict = {}
    names = list( seq_dict.keys() )

    query_pattern = r'OXX=[0-9]* *,[0-9]* *,[0-9]* *,[0-9]* *'

    for name in names:
        tax_ids = re.search( query_pattern, name ).group()
        ids     = tax_ids.split( '=' )[ 1 ]
        ids     = ids.split( ',' )
        for id in ids:
            if id not in id_dict:
                id_dict[ id ] = set()
            id_dict[ id ].add( name )
    return id_dict

def get_id_set( filename ):
    out_set = set()
    with open( filename, 'r' ) as open_file:
        for line in open_file:
            out_set.add( line.strip() )
    return out_set
            

if __name__ == '__main__':
    main()
