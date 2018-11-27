#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import sys                            # for exiting upon failure
from enum import Enum                 # for handling errors
import protein_oligo_library as oligo # for filling tax gaps

def main():

    parser = argparse.ArgumentParser( description = "Script to convert UniprotKB/Swiss-Prot "
                                                        "entries to a fasta file. The names for "
                                                        "each entry can be specified from the tags "
                                                        "found in each entry."
                                        )

    parser.add_argument( '-s', '--swiss', type = str,
                             help = "Name of file containing UniprotKB/Swiss-prot "
                                    "entries."
                           )
    parser.add_argument( '-o', '--output', type = str,
                             help = "Name of file to write outputs to, "
                                    "data in file will "
                                    "be formatted in the FASTA format."
                           )
    parser.add_argument( '-t', '--tags', action = 'append',
                             default = [ 'ID' ],
                             help = "Tags that will be included in each sequence name, "
                                    "note that id will always be collected and used for "
                                    "sequence names."
                                    "To include multiple tags, provide this argument multiple times, "
                                    "each time providing a tag to include. "
                                    "Note that if 'OC' is provided, then the "
                                    "taxdata flag must also be provided. "
                                    "Note that the OXX tag is not part of the "
                                    "uniprot standard, and both OXX/OX may be "
                                    "included. When OXX is included, each entry in the "
                                    "output will include the following information: "
                                    "OXX=original_id,species_id,genus_id,family_id. "
                                    "Note that if OXX is included, both the "
                                    "ranked_lineage and rank_map flags will need "
                                    "to be included. "
                                    "Possible values for tags include: "
                                    "AC, DE, DR, DT, FT, GN, ID (always included), "
                                    "KW, OC, OH, OS, OX, PE, RA, RC, RG, RL, RN, RP, RT, "
                                    "RX, SQ, OXX. "
                           )
    parser.add_argument( '-l', '--ranked_lineage', default = None, type = str,
                         help = "Map containing taxid|tax_info pairings "
                                "which can be used to create OC tags. "
                                "Note that is the OC tag is provided from "
                                "the command line, then this tag must be provided "
                                "as well. Inclusion of this argument without "
                                "the OC tag will be ignored."
                       )
    parser.add_argument( '-m', '--rank_map', default = None, type = str,
                         help = "Map containing taxid|tax_rank pairings, which will be "
                                "be used to annotate taxid taxonomic rank info. "
                                "Note that this argument is optional when combined with "
                                "the 'OC' tag. Note that this argument will be ignored "
                                "if 'OC' and '--ranked_lineage' are not also "
                                "provided. If provided, each OC tag will be "
                                "formatted 'OC=12345,FAMILY', where 'FAMILY' "
                                "is the taxonomic "
                                "rank identifier for the id."
                       )

    args = parser.parse_args()

    args_result = validate_args( args )

    if args_result != ArgResults.NO_ERR.value:
        report_error( args_result )
        sys.exit( 1 )
        
    tags = [ item.upper() for item in args.tags ]

    
    infile_parsers = Parser.create_parsers( [ ( get_taxdata_from_file, args.ranked_lineage ),
                                              ( get_rank_map_from_file, args.rank_map )
                                            ]
                                          )

    taxdata_from_cl, rank_map_from_cl = Parser.parse_list( infile_parsers )

    # parse the swisskb file
    db_parser = DBParser( args.swiss, args.output,
                          tags, taxdata = taxdata_from_cl,
                          rank_map = rank_map_from_cl
    )

    # convert the swisskb file to the FASTA format
    sequences = db_parser.parse()

    # write the output file to args.outfile
    write_outputs( args.output, sequences )
        

    sys.exit( 1 )


class DBParser:
    def __init__( self, db_filename, out_filename, tags_list, taxdata = None, rank_map = None ):
        self._tag_names = [ 'AC', 'DE', 'DR', 'DT',
                            'FT', 'GN', 'ID', 'KW',
                            'OC', 'OH', 'OS', 'OX',
                            'PE', 'RA', 'RC', 
                            'RL', 'RN', 'RP', 'RT',
                            'OXX'
                          ]
        self._db_filename  = db_filename
        self._out_filename = out_filename
        self._tags_list    = tags_list

        self._sequences    = list()
        self._taxdata      = taxdata
        self._rank_map     = rank_map

        if taxdata and rank_map:
            self._fixed_tax_data = oligo.fill_tax_gaps( taxdata, rank_map )
        else:
            self._fixed_tax_data = None

    def parse( self ):

        with open( self._db_filename, 'r' ) as open_file:
            seq_flag = False

            seqs = list()

            for line in open_file:
                split_line = line.split()
                tag_name = split_line[ 0 ]

                if tag_name != '//':
                    if seq_flag:
                        current_seq.add_seq_data( ''.join( split_line ) )
                    elif tag_name == 'ID':
                        current_seq = Sequence()
                        seqs.append( current_seq )

                        new_tag = DBParser.get_line_data(
                            tag_name,
                            split_line[ 1:: ],
                            self._tags_list,
                            self._taxdata,
                            self._rank_map
                        )

                        current_seq.add_tag( new_tag )

                    elif tag_name != 'SQ':
                        new_tag = DBParser.get_line_data(
                            tag_name,
                            split_line[ 1:: ],
                            self._tags_list,
                            self._taxdata,
                            self._rank_map
                        )

                        current_seq.add_tag( new_tag )

                        if tag_name == 'OX' and 'OXX' in self._tags_list:
                            new_tag = DBParser.get_line_data(
                                'OXX',
                                split_line[ 1:: ],
                                self._tags_list,
                                self._taxdata,
                                self._rank_map
                            )

                            current_seq.add_tag( new_tag )
                            
                    elif tag_name == 'SQ':
                        seq_flag = True

                else:
                    seq_flag = False

        return seqs
                            


    @staticmethod
    def get_line_data( str_tag_name, list_line, list_of_tags,
                       taxdata = None, rank_map = None, fixed_tax_data = None ):
        output_tag = None

        data_factory = TagDataFactory( taxdata, rank_map, fixed_tax_data )

        if str_tag_name in list_of_tags:
            new_tag = data_factory.create_tag( str_tag_name )
            new_tag.process( ' '.join( list_line ) )

            output_tag = new_tag

        return output_tag



class Sequence:
    def __init__( self ):
        self.tags     = {}
        self.seq_data = ''

    def add_tag( self, new_tag ):
        if new_tag:
            if new_tag.tag_type not in self.tags:
                self.tags[ new_tag.tag_type ] = list()
            self.tags[ new_tag.tag_type ].append( new_tag )

    def add_seq_data( self, data ):
        self.seq_data += data.strip()
            
    def __str__( self ):
        out_str = '>'
        out_str += str( self.tags[ 'ID' ][ 0 ] ).strip()

        for tag in self.tags:
            if 'ID' not in tag:
                for current in self.tags[ tag ]:
                    out_str += '%s' % str( current )
        out_str += "\n%s\n" % ( self.seq_data )

        return out_str

class TagDataFactory:
    SPACE      = ' '
    SEMICOLON  = ';'
    PERIOD     = '.'
    COMMA      = ','
    
    delimiters = {
                      'AC': SEMICOLON,
                      'DE': SEMICOLON,
                      'DR': PERIOD,
                      'DT': PERIOD,
                      'FT': SPACE,
                      'GN': SEMICOLON,
                      'ID': SPACE,
                      'KW': SEMICOLON,
                      'OC': SEMICOLON,
                      'OH': PERIOD,
                      'OS': PERIOD,
                      'OX': SEMICOLON,
                      'PE': SEMICOLON,
                      'RA': COMMA,
                      'RC': SEMICOLON,
                      'RG': '',
                      'RL': PERIOD,
                      'RN': '',
                      'RP': PERIOD,
                      'RT': SEMICOLON
                 }
    def __init__( self, taxdata, rank_map, fixed_tax_data = None ):
        self.taxdata        = taxdata
        self.rank_map       = rank_map
        self.fixed_tax_data = fixed_tax_data

    def create_tag( self, tag_name ):
        delimiters = TagDataFactory.delimiters

        if is_tax_tag( tag_name ):
            return_data = TaxTagData( tag_name, delimiters[ tag_name ] )
        elif tag_name == 'ID':
            return_data = IDTagData( tag_name, delimiters[ tag_name ] )
        elif tag_name == 'OC':
            return_data = TagData( tag_name, delimiters[ tag_name ] )
        elif tag_name == 'OXX':
            return_data = OXXTagData( tag_name, delimiters[ tag_name ], self.rank_map, self.fixed_tax_data )
        else:
            return_data = TagData( tag_name, delimiters[ tag_name ] )
        return return_data


def is_tax_tag( name ):
    return name == 'OH' or name == 'OX' 

class TagData:
    def __init__( self, tag_type, delimiter ):
        self.tag_type  = tag_type
        self.delimiter = delimiter
        self.data      = list()
        
    def process( self, line ):
        split_line = line.split( self.delimiter )

        for current_item in split_line:
            if len( current_item ) > 0:
                self.data.append( current_item.strip() )

    def __str__( self ):
        out_str = '' 
        for item in self.data:
            out_str += ' %s=%s' % ( self.tag_type, item )
        return out_str

class IDTagData( TagData ):
    def __init__( self, tag_name, delimiter ):
        super().__init__( tag_name, delimiter )

    def process( self, line ):
        self.data.append( line.split()[ 0 ] )

class TaxTagData( TagData ):
    def __init__( self, tag_name, delimiter ):
        super().__init__( tag_name, delimiter )

    def process( self, line ):
        tax_delimiter = 'NCBI_TaxID='
        split_line = line.split( tax_delimiter )

        if len( split_line[ 1 ] ) > 0:
            id_only = split_line[ 1 ].split( ';' )
            self.data.append( id_only[ 0 ] )

class OXXTagData( TagData ):
    reversed_table = None

    def __init__( self, tag_name, delimiter, taxdata, fixed_tax_data ):
        super().__init( tag_name, delimiter )
        self.taxdata = taxdata
        self.fixed_tax_data = fixed_tax_map

        if OXXTagData.reversed_table is None:
            OXXTagData.reversed_table = OXXTagData.reverse_table( taxdata )

    def process( self, line ):
        split_line = line.split( self.delimiter )

    @staticmethod
    def reverse_table( to_reverse ):
        out_table = {}
        for key, value in to_reverse.items():
            new_key   = value[ 0 ].strip()
            new_value = key.strip()
            out_table[ new_key ] = new_value 

        return out_table
        
class OCTagData( TagData ):
    def __init__( self, tag_name, delimiter, taxdata, tax_rank = None ):
        super().__init__( tag_name, delimiter )
        self.taxdata  = taxdata
        self.tax_rank = tax_rank

    def process( self, line ):
        split_line = line.split( self.delimiter )

        if len( split_line ) > 0: 
            for item in split_line:
                item = item.strip().lower()

                item = item.replace( '.', '' )

                try:
                    if self.tax_rank:
                        new_item = "%s,%s" % ( self.taxdata[ item ],
                                               self.tax_rank[ self.taxdata[ item ] ]
                                             )
                    else:
                        new_item = self.taxdata[ item ]
                    self.data.append( new_item )
                except KeyError:
                    pass # empty string

def write_outputs( outfile_name, seq_list ):
    with open( outfile_name, 'w' ) as out_file:
        for current_seq in seq_list:
            out_file.write( str( current_seq ) )

def validate_args( args_obj ):
    if 'OXX' in args_obj.tags and ( args_obj.ranked_lineage or \
                                    args_obj.rank_map is None ):
        return ArgResults.MISSING_TAXDATA_TAG.value
    return ArgResults.NO_ERR.value

def report_error( int_err_code ):
    err_codes = { 1: 'No Error',
                  2: 'OXX tag was provided by command '
                     'line, but ranked_lineage and rank_map '
                     'need to be provided.'
                }
    print( "ERROR: %s, program will exit..." % err_codes[ int_err_code ] )

class ArgResults( Enum ):
    BLANK                = 0,
    NO_ERR               = 1,
    MISSING_TAXDATA_TAG = 2

def get_taxdata_from_file( filename ):
    return_data = {}
    if filename:
        for key, value in oligo.get_taxdata_from_file( filename ).items():
            return_data[ str( key ) ] = value
        return return_data
    return None

class Parser:
    def __init__( self, parse_method, filename ):
        self.filename     = filename
        self.parse_method = parse_method

    def parse( self ):
        return self.parse_method( self.filename )

    @staticmethod
    def create_parser( parse_method, filename ):
        return Parser( parse_method, filename )

    @staticmethod
    def create_parsers( list_to_create ):
        out_parsers = list()
        for current_tuple in list_to_create:
            new_item = Parser.create_parser( current_tuple[ 0 ], current_tuple[ 1 ] )
            out_parsers.append( new_item )
        return out_parsers

    @staticmethod
    def parse_list( parser_list ):
        outputs_list = list()

        for item in parser_list:
            if item:
                outputs_list.append( item.parse() )
        return outputs_list

def process_line( str_line ):
    split_line = str_line.split( '|' )

    return split_line[ 1 ].strip().lower(), split_line[ 0 ].strip().lower()

def get_rank_map_from_file( filename ):
    rank_data = {}
    if filename:
        with open( filename, 'r' ) as open_file :
            for line in open_file:
                split_line = line.split( '|' )
                new_key = split_line[ 0 ].strip()
                new_val = split_line[ 1 ].strip()
                rank_data[ new_key ] = new_val 

        return rank_data
    return None

if __name__ == '__main__':
    main()
