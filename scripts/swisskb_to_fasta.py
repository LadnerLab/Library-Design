#!/usr/bin/env python3

import argparse                       # for parsing command-line arguments
import sys                            # for exiting upon failure
import protein_oligo_library as oligo # for file operations on sequences/fastas



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
                             default = [ 'id' ],
                             help = "Tags that will be included in each sequence name, "
                                    "note that id will always be collected and used for "
                                    "sequence names."
                                    "To include multiple tags, provide this argument multiple times, "
                                    "each time providing a tag to include. "
                                    "Possible values for tags include: "
                                    "AC, DE, DR, DT, FT, GN, ID (always included), "
                                    "KW, OC, OH, OS, OX, PE, RA, RC, RG, RL, RN, RP, RT, "
                                    "RX, SQ. "
                           )

    args = parser.parse_args()

    # parse the swisskb file
    db_parser = DBParser( args.swiss, args.output, args.tags )

    db_parser.parse()

    # convert the swisskb file to the FASTA format

    # write the output file to args.outfile
        

    sys.exit( 1 )


class DBParser:
    def __init__( self, db_filename, out_filename, tags_list ):
        self._tag_names = [ 'AC', 'DE', 'DR', 'DT',
                            'FT', 'GN', 'ID', 'KW',
                            'OC', 'OH', 'OS', 'OX',
                            'PE', 'RA', 'RC', 
                            'RL', 'RN', 'RP', 'RT'
                          ]
        self._tags = {}
                            
        self._db_filename  = db_filename
        self._out_filename = out_filename
        self._tags_list    = tags_list

    def parse( self ):

        tag_dict = self._tags

        with open( self._db_filename, 'r' ) as open_file:
            for line in open_file:
                split_line = line.split()
                tag_name = split_line[ 0 ]
                get_line_data( tag_name, split_line[ 1:: ], tags_list )

    @staticmethod
    def get_line_data( str_tag_name, list_line, list_of_tags ):
        output_list = list()

        data_factory = TagDataFactory()

        if str_tag_name in list_of_tags:
            new_tag = data_factory.create_tag( str_tag_name )
            new_tag.process( list_line )

            output_list.append( new_tag )


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
                      'OH': SEMICOLON,
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
    def __init__( self ):
        pass

    def create_tag( self, tag_name ):
        return_data = TagData( tag_name, delimiters[ tag_name ) )
        return return_data


if __name__ == '__main__':
    main()
