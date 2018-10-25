#!/usr/bin/env python3
import sys

def main():
    in_file = sys.argv[ 1 ]

    open_file = open( in_file )
    out_lines = list()

    counter = 0

    for line in open_file:
        line = line.split( '\t' )

        if counter == 0:
            out_lines.append( line )
            counter += 1
            continue

        if len( line[ 5 ] )  == 0:
            line[ 2 ] = str( int( line[ 2 ] ) - 1 )
        elif len( line[ 5 ].split( ',' ) ) > 0:
            line[ 2 ] = str( len( line[ 5 ].split( ',' ) ) )

        if len( line[ 6 ] )  == 0:
            line[ 3 ] = str( int( line[ 3 ] ) - 1 )
        elif len( line[ 6 ].split( ',' ) ) > 0:
            line[ 3 ] = str( len( line[ 6 ].split( ',' ) ) )

        if len( line[ 7 ] )  == 0:
            line[ 4 ] = str( int( line[ 4 ] ) - 1 )
        elif len( line[ 7 ].split( ',' ) ) > 0:
            line[ 4 ] = str( len( line[ 7 ].split( ',' ) ) )

        out_lines.append( line )

    out_file = open( "fixed_numbers.tsv", 'w' )
    for current in out_lines:
        out_file.write( '\t'.join( current )[:-1] + '\n' )
    out_file.close()



if __name__ == '__main__':
    main()
