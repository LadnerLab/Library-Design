#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = 'Plot coverage depth of the kmers for a reference in a design.' )
    arg_parser.add_argument( '-d', '--design', help = "Name of designed fasta file containing peptides." )
    arg_parser.add_argument( '-r', '--reference', help = "Name of reference fasta sequences." )
    arg_parser.add_argument( '-o', '--output', help = "Name of the directory to write plots to." )
    arg_parser.add_argument( '-k', '--kmer_size', help = "Size of kmer to use", type = int, default = 7 )

    args = arg_parser.parse_args()

    design_seqs = parse_fasta( args.design )
    ref_seqs    = parse_fasta( args.reference )

    do_kmers = lambda x: ( x.get_name(),
                           Collection( get_kmers_with_locs( x.get_seq(), args.kmer_size ) )
                         )
    
    ref_kmers = ref_seqs.apply(
                                 do_kmers,
                                 ref_seqs
                              ) 
    print( ref_kmers[ 0 ] )


def parse_fasta( fname ):
    names, sequences = oligo.read_fasta_lists( fname )
    sequences = [ Sequence( n, s ) for n, s in zip( names, sequences ) ] 

    return Collection( data = sequences ) 

def get_kmers_with_locs( str_seq, k ):
    xmers = set()

    start = 0
    end   = k

    while end <= len( str_seq ):
        xmer = str_seq[ start:end ]

        if 'X' not in xmer:
            xmers.add( Kmer( seq = xmer,
                             start = start + 1,
                             end = end
                           )
                     )

        start += 1
        end += k
    return xmers

    

class Collection:
        def __init__( self, data = None ):
            self.data = data if data != None else list()

        def __getitem__( self, idx ):
            return self.data[ idx ]

        def __setitem__( self, idx, value ):
            self.data[ idx ] = value

        def get_attr( self, attr ):
            return [ getattr( item ) for item in self.data ]

        def apply( self, func, what ):
            return [ func( item ) for item in what ]
            
        
class Sequence:
    def __init__( self, name = "", seq = "" ):
        self._name = name
        self._seq = seq

    def get_name( self ):
        return self._name

    def get_seq( self ):
        return self._seq

class Kmer( Sequence ):
    def __init__( self, name = "", seq = "", start = 0, end = 0 ):
        super().__init__( name = name, seq  = seq )
        self._start = start
        self._end = end

    def get_start( self ):
        return self._start

    def get_end( self ):
        return self._end

    def get_coordinates( self ):
        return ( self.get_start(), self.get_end() )

    def __hash__( self ):
        return hash( self._seq )

    def __eq__( self, other ):
        return self._seq == other.seq

    def __ne__( self, other ):
        return not( self.__eq__( other ) )

    def __lt__( self, other ):
        return self._start < other._start

if __name__ == '__main__':
    main()
