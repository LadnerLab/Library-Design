#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo
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

    args = arg_parser.parse_args()

    input_file = args.input

    output = args.output
    if not args.output:
        output = input_file + '_encodable.csv'
    line_width = args.line_size

    names, seqs = oligo.read_fasta_lists( input_file )

    oligo_len = len( seqs[ 0 ] )
    print( "Oligo length: ", oligo_len )
    allowed_length = line_width - oligo_len -1

    out_names = list()
    out_seqs  = list()

    for index, name in enumerate( names ):
        out_names.append( name[ 0 : allowed_length - 2].replace( ',', '~' )) # Commas not allowed in names
        out_seqs.append( seqs[ index ] )

    with open( output, 'w' ) as out_file:
        for index, name in enumerate( out_names ):
            out_file.write( "%s,%s\n" % ( name, seqs[ index ] ) )
        

if __name__ == '__main__':
    main()

