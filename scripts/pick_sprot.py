#!/usr/bin/env python3

import argparse                        # For parsing command-line arguments
import protein_oligo_library as oligo  # for operations on Fasta files


def main():
    arg_parser = argparse.ArgumentParser( description = "Parse representative output map to produce a sprot/trembl file" )

    arg_parser.add_argument( '-s', '--sprot', help = "Input sprot fasta to parse" )
    arg_parser.add_argument( '-t', '--trembl', help = "Input trembl file to parse" )
    arg_parser.add_argument( '-m', '--map_file', help = "Input map file to parse." )

    args = arg_parser.parse_args()

    out_sprot_name = args.map_file + "_sprot"
    out_trembl_name = args.map_file + "_trembl"




if __name__ == '__main__':
    main()
