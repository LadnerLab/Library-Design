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

    sys.exit( 1 )



if __name__ == '__main__':
    main()
