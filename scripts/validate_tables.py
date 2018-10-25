#!/usr/bin/env python3
import sys
import argparse

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = (
                                                            "Determines and outputs number of unique "
                                                            "sequences/species/genera/families in a fasta file"
                                                        )
                                        )

    arg_parser.add_argument( '-f',  '--fasta',    help = "Input fasta from which to gather data")
    arg_parser.add_argument( '-t', '--tax_db',   help = "Name of file containing mappings of taxids -> rank data")
    arg_parser.add_argument( '-g',  '--gap_file', help = "File containing mappings of taxid->rank for use in filling gaps" )

    arg_parser.add_argument( '--oligo_file', help = "Parsed oligo map file" )

    args = arg_parser.parse_args()

    taxid_dict = {}
    gap_dict   = {}

    missing_id_key = {
                       10969	: 444185,
                       11619	: 2169991,
                       11630	: 2169993,
                       11806	: 353765,
                       45218	: 2169996,
                       45222	: 2169994,
                       45709	: 2169992,
                       489502	: 10407,
                       587201	: 10255,
                       587202	: 10255,
                       587203	: 10255,
                       1173522	: 11723,
                       1554474	: 1511807,
                       1554482	: 1330068,
                       1554483	: 1330491,
                       1554492	: 1330066,
                       1554494	: 1307800,
                       1554498	: 1511784,
                       1559366	: 1513237,
                       1560037	: 1131483
                     }

    # parse and store tax_db and gap_file info
    with open( args.gap_file ) as gap_file:
        for line in gap_file:
            line = line.split( '|' )
            gap_dict[ line[ 0 ] ] = line[ 1 ].strip()

    with open( args.tax_db, 'r' ) as tax_db:
        for line in tax_db:
            line = line.split( '|' )
            taxid_dict[ line[ 0 ].strip() ] = [ item.strip() for item in line[ 1:: ] ]

    # Fill gaps
    for id, info in taxid_dict.items():
        if str(id) in gap_dict:
            if gap_dict[str(id)] == "SPECIES":
                taxid_dict[id][1] = taxid_dict[id][0]
                taxid_dict[id][0] = ""
            elif gap_dict[str(id)] == "GENUS":
                taxid_dict[id][2] = taxid_dict[id][0]
                taxid_dict[id][0] = ""
            elif gap_dict[str(id)] == "FAMILY":
                taxid_dict[id][3] = taxid_dict[id][0]
                taxid_dict[id][0] = ""

    num_seqs = 0
    unique_species  = set()
    unique_genera   = set()
    unique_families = set()


    names, sequences = oligo.read_fasta_lists( args.fasta )

    num_seqs = len( names )

    tax_ids = set( [ oligo.get_taxid_from_name( item ) \
                     for item in names 
                   ]
                 )

    missing_ids = set()

    for current in tax_ids:
        try:
            if int( current ) in missing_id_key:
                current = str( missing_id_key[ int( current ) ] )

            unique_species.add(   taxid_dict[ current ][ 1 ] )
            unique_genera.add(    taxid_dict[ current ][ 2 ] )
            unique_families.add(  taxid_dict[ current ][ 3 ] )
        except KeyError:
            missing_ids.add( current )
    

    unique_species  = [ item for item in unique_species if len( item ) > 0 ]
    unique_genera   = [ item for item in unique_genera if len( item ) > 0 ]
    unique_families = [ item for item in unique_families if len( item ) > 0 ]

    print( "Number of sequences:       %d" %  num_seqs )
    print( "Number of unique species:  %d" %  len( unique_species  ) )
    print( "Number of unique genera:   %d" %  len( unique_genera   ) )
    print( "Number of unique families: %d" %  len( unique_families ) )

    print( "Missing ids: %s" % ",".join( list( missing_ids ) ) )

    species_from_file  = set()
    genera_from_file   = set()
    families_from_file = set()

    with open( args.oligo_file, 'r' ) as oligo_file:
        counter = 0
        for line in oligo_file:
            if not counter:
                counter += 1
                continue 

            line = line.split( '\t' )

            for item in line[ 5 ].split( ',' ):
                if len( item.strip() ) > 0:
                    species_from_file.add( item.strip() )
            for item in line[ 6 ].split( ',' ):
                if len( item.strip() ) > 0:
                    genera_from_file.add( item.strip() )
            for item in line[ 7 ].split( ',' ):
                if len( item.strip() ) > 0:
                    families_from_file.add( item.strip() )

    print( "Number of species from oligo table:  %d" % len( species_from_file ) )
    print( "Number of genera from oligo table:   %d" % len( genera_from_file ) )
    print( "Number of families from oligo table: %d" % len( families_from_file ) )


    print( "Species missing from reference dataset:  %s "   % "|".join( species_from_file  - set( unique_species ) ) )
    print( "Genera missing from reference dataset:   %s "   % "|".join( genera_from_file   - set( unique_genera ) ) )
    print( "Families missing from refernece dataset: %s "   % "|".join( families_from_file - set( unique_families ) ) )

    print()

    print( "Species missing from oligo table:  %s "   % "|".join( set( unique_species )  - species_from_file ) )
    print( "Genera missing from oligo table:   %s "   % "|".join( set( unique_genera )   - genera_from_file ) )
    print( "Families missing from oligo table: %s "   % "|".join( set( unique_families ) - families_from_file ) )

    
if __name__ == '__main__':
    main()
