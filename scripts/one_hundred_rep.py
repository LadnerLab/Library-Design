#!/usr/bin/env python3
import argparse                        # For parsing command-line arguments
import protein_oligo_library as oligo  # for operations on Fasta files
def main():
    arg_parser = argparse.ArgumentParser( description = "Find and choose the 100% representative sequences from a fasta." )

    arg_parser.add_argument( '-f', '--fasta', help = "Fasta file to read input from" )
    arg_parser.add_argument( '-r', '--representatives', help = "Output file to write representatives to" )
    arg_parser.add_argument( '-m', '--map', help = "(Optional) Name of map to write the map of what "
                                                   "sequences were collapsed under what sequences."
                           )

    args = arg_parser.parse_args()

    # TODO validate args


    seqs_from_file = read_fasta( args.fasta )




if __name__ == '__main__':
    main()


class Sequence:
    def __init__( self, name, seq ):
        self.name = name
        self.seq  = seq

    def __len__( self ):
        return( len( self.seq ) )

    def __str__( self ):
        return ">%s\n%s\n" % ( self.name, self.seq )

    names, sequences = oligo.read_fasta_lists( args.fasta )

class SequenceFactory:
    def __init__( self ):
        pass

    def create_seq( self, name, sequence ):
        return Sequence( name, sequence )
    


            

