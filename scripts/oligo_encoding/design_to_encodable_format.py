#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo
def main():
    if len( sys.argv ) != 3:
        print( "USAGE: design_to_encodable_format.py in_file out_file")
        sys.exit( 1 )

    names, seqs = oligo.read_fasta_lists( sys.argv[ 1 ] )

    oligo_len = len( seqs[ 0 ] )
    print( "Oligo length: ", oligo_len )
    allowed_length = 128 - oligo_len -1

    out_names = list()
    out_seqs  = list()

    for index, name in enumerate( names ):
        out_names.append( name[0:allowed_length-2].replace( ',', '_' ))
        out_seqs.append( seqs[ index ])

    with open( sys.argv[ 2 ], 'w' ) as out_file:
        for index, name in enumerate( out_names ):
            out_file.write( "%s,%s\n" % (name, seqs[ index]))
        

if __name__ == '__main__':
    main()

