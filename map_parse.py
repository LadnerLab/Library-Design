#!/usr/bin/env python3
import argparse
import sys

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Script to parse information given by validate_design's epitope map output" )

    arg_parser.add_argument( '-m', '--map', help = "Name of file containing epitope map" )
    arg_parser.add_argument( '-o', '--output', default = "parsed_map.txt",
                             help = (
                                      "Name of file to write parsed output to."
                                      " Output will be written in a tab-delimited format"
                                    )
                           )
    arg_parser.add_argument( '-t', '--tax_db',
                             help = "Name of file containing mappings of taxonomic id -> rank data"
                           )
    arg_parser.add_argument( '-v', '--verbose',
                             help = "Flag to add if output should be written to STDOUT"
                           )

    args = arg_parser.parse_args()

    map_dict = None
    try:
        map_dict, seq_dict = parse_map( args.map )
    except ( IOError, OSError, TypeError ):
        print( "ERROR: An IO exception occurred when trying to "
               "open and parse map file"
             )
        sys.exit( 1 )
    except InputFormatFileError:
        print( "ERORR: The input file provided "
               "is not formatted correctly and cannot be "
               "parsed by this program."
             )
        sys.exit( 1 )

    try:
        tax_dict = oligo_to_tax( map_dict, args.tax_db )
    except ( IOError, OSError ):
        print( "ERROR: An IO exception occurred when trying "
               "to open and parse the taxonomic database file."
             )
        sys.exit( 1 )
    except InputFormatFileError:
        print( "ERROR: The taxonomic database file provided "
               "as input to the program is incorrectly formatted."
             )
        sys.exit( 1 )

    # By this point, our data can be safely assumed as valid,
    #so we don't have to do any more verification

    oligo_centric_table    = create_oligo_centric_table( tax_dict, map_dict )
    sequence_centric_table = create_sequence_centric_table( seq_dict )

    write_outputs( args.output, oligo_centric_table, sequence_centric_table )
    

def parse_map( file_name ):
    """
        Opens, reads, and parses a file containing an epitope map.
        
        On error, raises the appropriate error.
     
        On successful operation, returns a dictionary containing
        epitope: list of items mapping
    """
    READ_FLAG  = 'r'
    TAB_CHAR   = '\t'
    DELIMITER_CHAR = '~'

    oligo_dict = {}
    seq_dict   = {}
    open_file  = None

    split_line   = ""
    new_dict_key = ""
    new_dict_val = ""

    sequence_dict_key = ""
    sequence_dict_val = ""

    # try to open the input file
    try:
        open_file = open( file_name, READ_FLAG )
        # loop through each line in file
        for line in open_file:
            # split the line on the tab character
            try:

                split_line = line.split( TAB_CHAR )
                # Information for oligo-centric data
                new_dict_key = split_line[ 0 ]
                new_dict_val = split_line[ 1 ]

                #information for sequence-centric data
                sequence_dict_key = remove_loc_markers( split_line[ 0 ] )

                if sequence_dict_key in seq_dict.key():
                    seq_dict[ sequence_dict_key ] += 1
                else:
                    seq_dict[ sequence_dict_key ] = 0

            except ( IndexError, Exception ):
                raise InputFormatFileError

            if len( new_dict_key ) == 0 \
               or len( new_dict_val ) == 0:
                raise InputFormatFileError
            
            # String is formatted correctly 
            else:
                if new_dict_key in oligo_dict.keys():
                    oligo_dict[ new_dict_key ] = oligo_dict[ new_dict_key ] + \
                                                new_dict_val.strip().split( DELIMITER_CHAR )
                else:
                    oligo_dict[ new_dict_key ] = new_dict_val.strip().split( DELIMITER_CHAR )
                    
    except:
        raise

    open_file.close()

    return oligo_dict, seq_dict

def oligo_to_tax( input_dict, tax_data_file ):
    """
       Parses data in input dict, and returns a dictionary
       containing seq_name: [ tax ranks ] for each entry
       in the dictionary
    """
    output_dict = {}
    taxid_dict = oligo.get_taxdata_from_file( tax_data_file )

    dict_keys = input_dict.keys()

    output_dict[ 'NoID' ] = list()

    for item in dict_keys:
        output_dict[ item ] = list()

        for current in input_dict[ item ]:
            current_tax_id = oligo.get_taxid_from_name( current ) 

            if current_tax_id \
                 and int( current_tax_id ) in taxid_dict.keys():

                current_tax_id = int( current_tax_id )
                output_dict[ item ].append( taxid_dict[ current_tax_id ] )

            else:
                output_dict[ 'NoID' ].append( current )

    return output_dict

def remove_loc_markers( input_str ):
    """
        Removes the location markers from a name in the epitope map file

        Input is the name of a sequences followed by int_int
     
        e.g. Seq1_23_45
    """
    UNDERSCORE_CHAR = '_'

    split_str = input_str.split( UNDERSCORE_CHAR )

    if len( split_str ) < 3:
        raise InputFormatFileError

    split_str.pop()
    split_str.pop()

    return '_'.join( split_str )

class InputFormatFileError( Exception ):
    pass

def create_oligo_centric_table( tax_dict, map_dict ):

    out_str = (  "Oligo Name\tNum Sequences Share 7-mer\tNum Species Share 7-mer\t"
                 "Num Genera Share 7-mer\t"
                 "Num Families Share 7-mer\n"
              )
    oligo_names     = map_dict.keys()
    num_oligos      = len( oligo_names )
    species_shared  = set()

    species_total   = 0
    genus_total     = 0
    family_total    = 0
    sequences_total = 0

    for current_oligo in oligo_names:
        current_species = get_num_items_at_rank( tax_dict[ current_oligo ],
                                                 oligo.Rank.SPECIES.value
                                               )
        current_genus   = get_num_items_at_rank( tax_dict[ current_oligo ],
                                                 oligo.Rank.GENUS.value
                                               )

        current_family  =  get_num_items_at_rank( tax_dict[ current_oligo ],
                                                  oligo.Rank.FAMILY.value
                                                )
        current_entry = current_oligo
        current_entry += "\t%d\t" % len( map_dict[ current_oligo ] )
        current_entry += "%d\t"   %  current_species
        current_entry += "%d\t"   %  current_genus
        current_entry += "%d\t"   %  current_family

        sequences_total += len( map_dict[ current_oligo ] )
        species_total   += current_species
        genus_total     += current_genus
        family_total    += current_family

        out_str         += current_entry + "\n"

    out_str += "Average\t%.2f\t%.2f\t%.2f\t%.2f\n" % ( ( sequences_total / num_oligos ),
                                                 ( species_total / num_oligos ),
                                                 ( genus_total   / num_oligos ),
                                                 ( family_total  / num_oligos )
                                                     )

    return out_str 

def get_num_items_at_rank( tax_list, rank ):
    shared_items = set()

    for current in tax_list:
        if len( current[ rank ] ) > 0:
            shared_items.add( current[ rank ] )
    return len( shared_items )

def write_outputs( out_file, oligo_centric, sequence_centric ):
    WRITE_FLAG = "w"
    EXTENSION  = ".tsv"

    oligo_file = open( out_file + "oligo_table" + EXTENSION,
                       WRITE_FLAG
                     )
    oligo_file.write( oligo_centric )
    oligo_file.close()

    sequence_file = open( out_file + "sequence_table" + EXTENSION,
                          WRITE_FLAG
                        )
    sequence_file.write( sequence_centric )
    sequence_file.close()   
    
if __name__ == '__main__':
    main()
