#!/usr/bin/env python3
import argparse                        # For parsing command-line arguments
import protein_oligo_library as oligo  # for operations on Fasta files
def main():
    arg_parser = argparse.ArgumentParser( description = "Find and choose the 100% representative sequences from a fasta." )

    arg_parser.add_argument( '-f', '--fasta', help = "Fasta file to read input from" )
    arg_parser.add_argument( '-r', '--representatives', help = "Output file to write representatives to" )
    arg_parser.add_argument( '-m', '--map_file', help = "(Optional) Name of map to write the map of what "
                                                   "sequences were collapsed under what sequences."
                           )

    args = arg_parser.parse_args()

    # TODO validate args


    input_parser = FastaParser( args.fasta )
    input_seqs   = input_parser.parse()
    indexer      = SortIndexer( len )

    # get the 100% reps for each sequence
    final_seqs, map_out = get_one_hundred_reps( input_seqs, indexer, args.map_file != None )

    print( "Number of seqs in original: %d" % len( input_seqs ) )
    print( "Number of seqs in output:   %d" % len( final_seqs ) )

    # write the output fasta and map file out
    write_outputs( args.representatives, final_seqs, args.map_file, map_out )

class Sequence:
    def __init__( self, name, seq ):
        self.name = name
        self.seq  = seq

    def __len__( self ):
        return( len( self.seq ) )

    def __eq__( self, other ):
        return self.seq == other.seq

    def __hash__( self ):
        return hash( self.seq )

    def __str__( self ):
        return ">%s\n%s\n" % ( self.name, self.seq )


class SequenceFactory:
    def __init__( self ):
        pass

    def create_seq( self, name, sequence ):
        return Sequence( name, sequence )
    def create_seq_list( self, names, sequences ):
        out_list = list()

        for index in range( len( names ) ):
            current_name     = names[ index ]
            current_sequence = sequences[ index ]

            out_list.append( self.create_seq( current_name, current_sequence ) )
        return out_list
    

class FastaParser:
    def __init__( self, filename ):
        self.filename    = filename
        self.seq_factory = SequenceFactory()

    def parse( self ):

        names, sequences = oligo.read_fasta_lists( self.filename )

        return self.seq_factory.create_seq_list( names, sequences )
            
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

def get_one_hundred_reps( seq_list, indexer, do_map = False ):
    seqs_set     = set()
    out_seqs     = list()

    out_map      = None

    if do_map:
        out_map = {}

    unique_seqs = get_unique_sequences( seq_list )

    indexed_seqs = indexer.index( unique_seqs, reverse = False )
    seqs_set     = set( indexed_seqs )
    seqs_str_set = set( [ item.seq for item in indexed_seqs ] )
    out_seqs     = set()

    index = 0
    combined_string = ''.join( seqs_str_set )

    for current_seq in indexed_seqs:
        seqs_str_set.remove( current_seq.seq )
        combined_string = ''.join( seqs_str_set )

        if current_seq.seq not in combined_string:
            out_seqs.add( current_seq )


    return list( out_seqs ), out_map

def get_unique_sequences( seq_list ):
    seq_fact = SequenceFactory()

    names_list = [ item.name for item in seq_list ]
    seq_list   = [ item.seq for item in seq_list ]

    new_names, new_seqs = oligo.get_unique_sequences( names_list, seq_list )

    return seq_fact.create_seq_list( new_names, new_seqs )

def write_outputs( seq_file, seq_list, map_file, out_map ):
    if seq_file:
        with open( seq_file, 'w' ) as open_file:
            for current_seq in seq_list:
                open_file.write( str( current_seq ) )
    if map_file:
        with open( map_file, 'w' ) as open_file:
            for seq, values in out_map.items():
                out_str = "%s\t%s\n" % ( seq.name, ','.join( [ item.name for item in values ] ) )
                open_file.write( out_str )
        

if __name__ == '__main__':
    main()


