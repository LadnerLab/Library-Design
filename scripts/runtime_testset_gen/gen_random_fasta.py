#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse
import random


def main():
    argp = argparse.ArgumentParser( description = "Randomly generate a fasta "
                                                  "containing AA sequences."
                                  )
    argp.add_argument( '-f', '--filename' )
    argp.add_argument( '-l', '--length', type = int )
    argp.add_argument( '-g', '--generate',
                       type = int
                     )
    argp.add_argument( '-p', '--prefix', help = "The prefix for each "
                                                "generated sequence. ",
                       default = "sequence"
                     )

    args = argp.parse_args()


if __name__ == '__main__':
    main()
