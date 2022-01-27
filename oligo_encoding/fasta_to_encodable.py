#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo
import math
import argparse
def main():
    arg_parser = argparse.ArgumentParser( "Transform a fasta file into a file of the format NAME,SEQ, note that names will be "
                                          "trimmed to fit line width requirements. Note that any commas present in the input names "
                                          "will become the '~' character upon output."
                                        )
    arg_parser.add_argument( '-i', '--input', help = "Name of file to transform, file must be in FASTA format." )
    arg_parser.add_argument( '-o', '--output', help = "Name of file to write output to, if not provided defaults to"
                                                      "input_encodable, where input is the name of input file provided."
                           )
    arg_parser.add_argument( '-s', '--line_size', help = "Max line width for output, names will be trimmed to fit these constraints. [128]",
                             type = int, default = 128
                           )
    arg_parser.add_argument( '-p', '--prefix', help = "Prefix for names that will be written to the encoded format. " )
    arg_parser.add_argument( '-m', '--map_file', help = "File to write tab-delimited pairings of original name to prefix name and number." )

    args = arg_parser.parse_args()

    input_file = args.input

    output = args.output
    if not args.output:
        output = input_file + '_encodable.csv'
    line_width = args.line_size

    original_names_and_seqs = fasta_to_dict( input_file )
    prefixed_names          = gen_prefixed_names( original_names_and_seqs.keys(), args.prefix )

    write_map( args.map_file, prefixed_names )

    prefixed_names_and_seqs = set_prefixed_names_to_seqs( original_names_and_seqs, prefixed_names )
    write_output( args.output, prefixed_names_and_seqs, args.line_size )

def write_output( filename, names_dict, line_size ):
    with open( filename, 'w' ) as open_file:
        for name, seq in names_dict.items():
            if len( name ) + len( seq ) + 2 >= line_size:
                print( "WARNING: %s,%s is too long of a line!" % ( name, seq ) )
            open_file.write( "%s,%s\n" % ( name, seq ) )
    
def set_prefixed_names_to_seqs( original_seq_dict, prefixed_name_dict ):
    out_dict = {}
    for name, seq in original_seq_dict.items():
        out_dict[ prefixed_name_dict[ name ] ] = seq
    return out_dict
    
def write_map( filename, names_dict ):
    with open( filename, 'w' ) as out_file:
        for original_name, prefixed_name in names_dict.items():
            out_file.write( "%s\t%s\n" % ( original_name, prefixed_name ) )

def fasta_to_dict( filename ):
    out_dict = {}
    names, sequences = oligo.read_fasta_lists( filename )
    for index, name in enumerate( names ):
        out_dict[ name ] = sequences[ index ]
    return out_dict

def gen_prefixed_names( names_list, prefix ):
    out_names = {}
    num_digits = int( math.log10( len( names_list ) ) + 1 )
    for index, current_name in enumerate( names_list ):
        digit_str = str( index ).zfill( num_digits )
        prefixed_name = "%s_%s" % ( prefix, digit_str )
        out_names[ current_name ] = prefixed_name
    return out_names

if __name__ == '__main__':
    main()

