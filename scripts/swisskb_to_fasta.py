#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import sys                            # for exiting upon failure
from enum import Enum                 # for handling errors



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
                                    "Possible values for tags include: "
                                    "AC, DE, DR, DT, FT, GN, ID (always included), "
                                    "KW, OC, OH, OS, OX, PE, RA, RC, RG, RL, RN, RP, RT, "
                                    "RX, SQ. "
                           )
    parser.add_argument( '-l', '--ranked_lineage', default = None, type = str,
                         help = "Map containing taxid|tax_info pairings "
                                "which can be used to create OC tags. "
                                "Note that is the OC tag is provided from "
                                "the command line, then this tag must be provided "
                                "as well. Inclusion of this argument without "
                                "the OC tag will be ignored."
                       )

    args = parser.parse_args()

    args_result = validate_args( args )

    if args_result != ArgResults.NO_ERR.value:
        report_error( args_result )
        sys.exit( 1 )
        
    tags = [ item.upper() for item in args.tags ]

    taxdata_from_cl = get_taxdata_from_file( args.ranked_lineage )

    # parse the swisskb file
    db_parser = DBParser( args.swiss, args.output, tags, taxdata = taxdata_from_cl )

    # convert the swisskb file to the FASTA format
    sequences = db_parser.parse()

    # write the output file to args.outfile
    write_outputs( args.output, sequences )
        

    sys.exit( 1 )


class DBParser:
    def __init__( self, db_filename, out_filename, tags_list, taxdata = None ):
        self._tag_names = [ 'AC', 'DE', 'DR', 'DT',
                            'FT', 'GN', 'ID', 'KW',
                            'OC', 'OH', 'OS', 'OX',
                            'PE', 'RA', 'RC', 
                            'RL', 'RN', 'RP', 'RT'
                          ]
        self._db_filename  = db_filename
        self._out_filename = out_filename
        self._tags_list    = tags_list

        self._sequences    = list()
        self._taxdata     = taxdata

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
                            self._taxdata
                        )

                        current_seq.add_tag( new_tag )

                    elif tag_name != 'SQ':
                        new_tag = DBParser.get_line_data(
                            tag_name,
                            split_line[ 1:: ],
                            self._tags_list,
                            self._taxdata
                        )

                        current_seq.add_tag( new_tag )
                    elif tag_name == 'SQ':
                        seq_flag = True

                else:
                    seq_flag = False

        return seqs
                            


    @staticmethod
    def get_line_data( str_tag_name, list_line, list_of_tags, taxdata = None ):
        output_tag = None

        data_factory = TagDataFactory( taxdata )

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
    def __init__( self, taxdata ):
        self.taxdata = taxdata

    def create_tag( self, tag_name ):
        if is_tax_tag( tag_name ):
            return_data = TaxTagData( tag_name, TagDataFactory.delimiters[ tag_name ] )
        elif tag_name == 'ID':
            return_data = IDTagData( tag_name, TagDataFactory.delimiters[ tag_name ] )
        elif tag_name == 'OC':
            return_data = OCTagData( tag_name, TagDataFactory.delimiters[ tag_name ], self.taxdata )
        else:
            return_data = TagData( tag_name, TagDataFactory.delimiters[ tag_name ] )
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

class OCTagData( TagData ):
    def __init__( self, tag_name, delimiter, taxdata ):
        super().__init__( tag_name, delimiter )
        self.taxdata = taxdata

    def process( self, line ):
        split_line = line.split( self.delimiter )

        if len( split_line ) > 0: 
            for item in split_line:
                item = item.strip().lower()
                if '.' in item:
                    item = item.replace( '.', '' )

                try:
                    self.data.append( self.taxdata[ item ] )
                except KeyError:
                    pass # empty string

def write_outputs( outfile_name, seq_list ):
    with open( outfile_name, 'w' ) as out_file:
        for current_seq in seq_list:
            out_file.write( str( current_seq ) )

def validate_args( args_obj ):
    if 'OC' in args_obj.tags and args_obj.ranked_lineage is None:
        return ArgResults.MISSING_TAXDATA_TAG.value
    return ArgResults.NO_ERR.value

def report_error( int_err_code ):
    err_codes = { 1: 'No Error',
                  2: 'OC tag was provided by command '
                     'line, but no taxdata file was provided'
                }
    print( "ERROR: %s, program will exit..." % err_codes[ int_err_code ] )

class ArgResults( Enum ):
    BLANK                = 0,
    NO_ERR               = 1,
    MISSING_TAXDATA_TAG = 2

def get_taxdata_from_file( filename ):
    return_data = {}
    if filename:
        with open( filename, 'r' ) as ranked_lineage:
            for current_line in ranked_lineage:
                key, value = process_line( current_line )
                return_data[ key ] = value

        return return_data
    return None

def process_line( str_line ):
    split_line = str_line.split( '|' )

    return split_line[ 1 ].strip().lower(), split_line[ 0 ].strip().lower()
    

if __name__ == '__main__':
    main()
