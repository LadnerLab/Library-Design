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


    input_parser = FastaParser( args.fasta )
    input_seqs   = input_parser.parse()
    indexer      = SortIndexer( len )

    # get the 100% reps for each sequence
    final_seqs, map_out = get_one_hundred_reps( input_seqs, indexer, args.map )

    # write the output fasta and map file out

class Sequence:
    def __init__( self, name, seq ):
        self.name = name
        self.seq  = seq

    def __len__( self ):
        return( len( self.seq ) )

    def __str__( self ):
        return ">%s\n%s\n" % ( self.name, self.seq )


class SequenceFactory:
    def __init__( self ):
        pass

    def create_seq( self, name, sequence ):
        return Sequence( name, sequence )
    

class FastaParser:
    def __init__( self, filename ):
        self.filename    = filename
        self.seq_factory = SequenceFactory()

    def parse( self ):
        out_seqs = list()

        names, sequences = oligo.read_fasta_lists( self.filename )

        for index in range( len( names ) ):
            new_name = names[ index ]
            new_seq  = sequences[ index ]
            out_seqs.append( self.seq_factory.create_seq( new_name, new_seq ) )

        return out_seqs
            
class Indexer:
    def __init__( self ):
        pass
    def index( self ):
        pass

class SortIndexer( Indexer ):
    def __init__( self, sort_key = len ):
        self.sort_key = sort_key

    def index( self, in_list, reverse = False ):
        return( sorted( in_list, key = self.sort_key, reverse = reverse ) )

    
if __name__ == '__main__':
    main()


