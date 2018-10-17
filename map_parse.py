#!/usr/bin/env python3
import argparse
import sys

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Script to parse information given by validate_design's epitope map output" )

    arg_parser.add_argument( '-m', '--map', help = "Name of file containing epitope map" )
    arg_parser.add_argument( '-o', '--output', default = "parsed_map",
                             help = (
                                      "Identifier to preface output files."
                                      " Output will be written in a tab-delimited format"
                                    )
                           )
    arg_parser.add_argument( '-t', '--tax_db',
                             help = "Name of file containing mappings of taxonomic id -> rank data"
                           )
    arg_parser.add_argument( '-v', '--verbose',
                             help = "Flag to add if output should be written to STDOUT"
                           )
    arg_parser.add_argument( '-g', '--gap_file',
                             help = ( "File containing mappings of taxid->rank for "
                                      "use in filling gaps."
                                    )
                           )
    arg_parser.add_argument( '-r', '--reference',
                             help = "Fasta reference dataset used to create library"
                           )

    args = arg_parser.parse_args()

    oligo_centric_table    = None
    sequence_centric_table = None
    species_centric_table  = None

    map_dict = None
    gap_dict = parse_gaps( args.gap_file )

    try:
        map_dict, epitope_dict, oligo_seq_dict = parse_map( args.map )
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
        tax_dict, taxid_dict = oligo_to_tax( map_dict, args.tax_db )
        sequence_dict        = oligo.sequence_dict_from_file( args.reference )
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

    # species_centric_table  = create_species_centric_table( map_dict, taxid_dict, oligo_seq_dict, tax_dict, gap_dict )
    # oligo_centric_table    = create_oligo_centric_table( tax_dict, map_dict, taxid_dict, gap_dict )
    # sequence_centric_table = create_sequence_centric_table_from_oligos( epitope_dict, oligo_seq_dict, gap_dict )

    create_sequence_centric_table( sequence_dict, epitope_dict, gap_dict )

    write_outputs( args.output, oligo_centric_table, sequence_centric_table )
    

def parse_map( file_name ):
    """
        Opens, reads, and parses a file containing an epitope map.
        
        On error, raises the appropriate error.
     
        On successful operation, returns a dictionary containing
        epitope: list of items mapping
        
        On successful operation, returns a dictionary containing
        sequence: number of occurences mapping
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

    oligo_list     = list()
    seq_list       = set()
    oligo_seq_dict = {}

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

                if sequence_dict_key in oligo_seq_dict.keys():
                    oligo_seq_dict[ sequence_dict_key ] += 1
                else:
                    oligo_seq_dict[ sequence_dict_key ] = 1

                for current_item in split_line[ 1 ].split( DELIMITER_CHAR ):
                    if current_item in oligo_seq_dict.keys():
                        oligo_seq_dict[ current_item ] += 1
                    else:
                        oligo_seq_dict[ current_item ] = 1

                if sequence_dict_key in seq_dict.keys():
                    current_entry = seq_dict[ sequence_dict_key ]
                    current_entry[ 0 ] += 1
                    current_entry[ 1 ] += len( split_line[ 1 ].split( DELIMITER_CHAR ) )
                else:
                    seq_dict[ sequence_dict_key ] = [ 1, len( split_line[ 1 ].split( DELIMITER_CHAR ) ) ]

            except ( IndexError ):
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

    return oligo_dict, seq_dict, oligo_seq_dict

def oligo_to_tax( input_dict, tax_data_file, gap_dict = None ):
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

    return output_dict, taxid_dict

def remove_loc_markers( input_str ):
    """
        Removes the location markers from a name in the epitope map file

        Input is the name of a sequences followed by int_int
     
        e.g. Seq1_23_45 -> Seq1
    """
    UNDERSCORE_CHAR = '_'

    split_str = input_str.split( UNDERSCORE_CHAR )

    if len( split_str ) < 3:
        raise InputFormatFileError

    split_str.pop()
    split_str.pop()

    return UNDERSCORE_CHAR.join( split_str )

class InputFormatFileError( Exception ):
    pass

def create_oligo_centric_table( tax_dict, map_dict, taxid_dict, gap_dict = None ):

    out_str = (  "Oligo Name\tNum Sequences Share 7-mer\tNum Species Share 7-mer\t"
                 "Num Genera Share 7-mer\t"
                 "Num Families Share 7-mer\t"
                 "Species Covered by 7-mer\t"
                 "Genera Covered by 7-mer\t"
                 "Families Covered by 7-mer\n"
              )
    oligo_names     = map_dict.keys()
    num_oligos      = len( oligo_names )
    species_shared  = set()

    species_total   = set()
    genus_total     = set()
    family_total    = set()
    sequences_total = 0

    for current_oligo in oligo_names:
        current_taxid   = oligo.get_taxid_from_name( current_oligo )
        current_entry = "%s\t%s\t" % ( current_oligo, len( tax_dict[ current_oligo ] ) )

        try:
            if current_taxid:
                rank_data = gap_dict[ current_taxid ].strip()

                rank_val = oligo.Rank[ rank_data ].value

            current_species, \
            current_genus, \
            current_family = get_items_at_rank_from_seqs( map_dict[ current_oligo ],
                                                          taxid_dict, 
                                                          gap_dict,
                                                          oligo.Rank.FAMILY.value
                                       )  

            current_entry += "%d\t%d\t%d\t" % ( len( current_species ),
                                                len( current_genus ),
                                                len( current_family )
                                              )

            current_entry += "%s\t" % ",".join( [ item for item in current_species \
                                                if len( item ) > 0 ]
                                              )
            current_entry += "%s\t" % ",".join( [ item for item in current_genus \
                                                  if len( item ) > 0 ]
                                              )
            current_entry += "%s" % ",".join( [ item for item in current_family \
                                                  if len( item ) > 0 ]
                                              )

            species_total   |= set( current_species )
            genus_total     |= set( current_genus )
            family_total    |= set( current_family )

            
            sequences_total += len( map_dict[ current_oligo ] )
            out_str         += current_entry + "\n"

        except KeyError:
            pass

    return out_str 

def get_items_at_rank( tax_list, rank ):
    shared_items = set()
    num_items    = 0

    if not tax_list:
        num_items = 1

    for current in tax_list:
        if len( current[ rank ] ) > 0:
            shared_items.add( current[ rank ] )
        else:
            num_items += 1
    return ( shared_items, len( shared_items ) + num_items ) # +1 because each rank is represented by at least this sequence

def write_outputs( out_file, oligo_centric, sequence_centric ):
    WRITE_FLAG = "w"
    EXTENSION  = ".tsv"

    if oligo_centric:
        oligo_file = open( out_file + "_oligo_table" + EXTENSION,
                           WRITE_FLAG
                         )
        oligo_file.write( oligo_centric )
        oligo_file.close()

    if sequence_centric:
        sequence_file = open( out_file + "_sequence_table" + EXTENSION,
                              WRITE_FLAG
                            )
        sequence_file.write( sequence_centric )
        sequence_file.close()   
    
def create_sequence_centric_table_from_oligos( seq_dict, oligo_seq_dict, gap_dict = None ):
    dict_keys  = seq_dict.keys()
    out_string = ( "Sequence Name\t"
                   "Number Oligos Seq Contrib. to Design\t"
                   "Number Seqs share 7-mer\t"
                   "Number Oligos Share 7-mer\n"
                 )

    for item in dict_keys:
        out_string += "%s\t%d\t%d\t%d\n" % ( item, seq_dict[ item ][ 0 ],
                                             seq_dict[ item ][ 1 ],
                                             oligo_seq_dict[ item ]
                                           )
        
    return out_string

def create_species_centric_table( map_dict, taxid_dict, oligo_seq_dict, tax_dict, gap_dict = None ):
    dict_keys      = map_dict.keys()
    seq_tax_dict   = {}
    oligo_dict     = {}

    out_string = ( "Species Name\t"
                   "Number Oligos Species Contrib. to Design\t"
                   "Number Species share 7-mer\t"
                   "Number Oligo Species Share 7-mer\t"
                   "Species that Share 7-mer\n"
                 )

    for item in dict_keys:
        current_taxid = oligo.get_taxid_from_name( item )

        if current_taxid:
            try:
                rank_data     = gap_dict[ current_taxid ].strip()
                rank_val      = oligo.Rank[ rank_data ].value

                if rank_val >= oligo.Rank.FAMILY.value: 

                    species = taxid_dict[ int( current_taxid ) ][ oligo.Rank.SPECIES.value ]
                                                       
                elif rank_val == oligo.Rank.SPECIES.value:

                    species = taxid_dict[ int( current_taxid ) ][ 0 ]
                elif rank_val == oligo.Rank.GENUS.value:
                    print( rank_val )
                    
            except KeyError:
                pass

            species = species.strip() 

            if species not in seq_tax_dict.keys():
                seq_tax_dict[ species ] = [ 1, set() ]
            else:
                seq_tax_dict[ species ][ 0 ] += 1

            for current_value in map_dict[ item ]:

                current_taxid = oligo.get_taxid_from_name( current_value )
                if current_taxid:
                    try:
                        rank_data     = gap_dict[ current_taxid ].strip()
                        rank_val      = oligo.Rank[ rank_data ].value

                        if rank_val >= oligo.Rank.FAMILY.value: 

                            current_species = taxid_dict[ int( current_taxid ) ][ oligo.Rank.SPECIES.value ]
                                                               
                        elif rank_val == oligo.Rank.SPECIES.value:

                            current_species = taxid_dict[ int( current_taxid ) ][ 0 ]
                        elif rank_val == oligo.Rank.GENUS.value:
                            print( rank_val )
                            
                    except KeyError:
                        pass


                    seq_tax_dict[ species ][ 1 ].add( current_species )

    for item in seq_tax_dict.keys():
        for current_species in seq_tax_dict[ item ][ 1 ]:
            if current_species not in oligo_dict.keys():
                oligo_dict[ current_species ] = 0

            oligo_dict[ current_species ] += 1

    for item in seq_tax_dict.keys():
        if len( item ) > 0:
            out_string += "%s\t%d\t%d\t%d\t%s\n" % ( item, seq_tax_dict[ item ][ 0 ],
                                                 len( seq_tax_dict[ item ][ 1 ] ),
                                                 oligo_dict[ item ],
                                                 ",".join( [ species for species \
                                                             in seq_tax_dict[ item ][ 1 ] \
                                                             if len( species ) > 0
                                                           ]
                                                         )
                                                   )

    print( out_string )
    return out_string

        
def parse_gaps( gap_file ):
    return_dict = {}

    return_val = None

    open_file = None

    if gap_file:
        open_file = open( gap_file, 'r' )

        for line in open_file:
            split_line = line.split( '|' )
            return_dict[ split_line[ 0 ].strip() ] = split_line[ 1 ]

        open_file.close() 

        return_val = return_dict
        
    return return_val 

def get_items_at_rank_from_seqs( oligo_items,
                                 taxid_dict,
                                 gap_dict, 
                                 oligo_rank
                               ):

    out_species  = set()
    out_genera   = set()
    out_families = set()

    for current_name in oligo_items:

        current_id = oligo.get_taxid_from_name( current_name ).strip()
        
        try:
            if oligo_rank >= oligo.Rank.FAMILY.value:
                out_families.add( taxid_dict[ int( current_id ) ][ oligo_rank ] )
                out_genera.add( taxid_dict[ int( current_id ) ][ oligo.Rank.GENUS.value ] )
                out_species.add( taxid_dict[ int( current_id ) ][ oligo.Rank.SPECIES.value ] )

            elif oligo_rank == oligo.Rank.GENUS.value:
                out_genera.add( taxid_dict[ int( current_id ) ][ 0 ] )

            elif oligo_rank == oligo.Rank.SPECIES.value:
                out_species.add( taxid_dict[ int( current_id ) ][ 0 ] )

        except KeyError:
            pass
            
    return out_species, out_genera, out_families

def create_sequence_centric_table( ref_dict, epitope_dict, gap_dict ):
    out_string = ( "Sequence Name\t"
                   "Oligos that Share 7-mer with Sequence\t"
                   "Oligos Share 7-mer with only sequence\t"
                 )

    oligo_share_dict = {}

    for current_seq in ref_dict.keys():
        for current_oligo in epitope_dict:
            if current_seq in epitope_dict[ current_oligo ]:
                if current_seq not in oligo_share_dict.keys():
                    oligo_share_dict[ current_seq ] = 1
                oligo_share_dict[ current_seq ] += 1

        out_string += current_seq
    return out_string

if __name__ == '__main__':
    main()
