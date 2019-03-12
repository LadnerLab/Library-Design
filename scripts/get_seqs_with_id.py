#!/usr/bin/env python3
import argparse
import sys


def main():
    argparser = argparse.ArgumentParser( description = "Retrieve all sequences in a file that "
                                                       "have contain a specified taxonomic id."
                                       )
    argparser.add_argument( '-i', '--ids', help = "Name of file containing one "
                                                  "id per line"
                          )
    argparser.add_argument( '-f', '--fasta', help = "Name of fasta file whose sequences "
                                                    "will be searched for the specified ids."
                          )
    argparser.add_argument( '-t', '--tags', help = "Comma-separated list of tags that specify "
                                                   "taxonomic id filed in sequence names. [OXX]",
                            default = "OXX"
                          )

    args = argparser.parse_args()
                                                          
if __name__ == '__main__':
    main()
