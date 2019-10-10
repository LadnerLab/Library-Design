#!/usr/bin/env python3
import argparse
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
import protein_oligo_library as oligo
import os
import numpy as np

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
                           get_kmers_with_locs( x.get_seq(), args.kmer_size )
                         )
    
    ref_kmers = ref_seqs.apply(
                                 do_kmers,
                                 ref_seqs
                              ) 
    des_kmers = kmer_dict_wcounts( design_seqs, args.kmer_size )

    positional_counts = list()
    for ref_seq in ref_kmers:
        name, kmers = ref_seq

        tup = ( name, get_positional_counts( kmers, des_kmers ) ) 
        positional_counts.append( tup )

    create_figures( positional_counts, args.output, args.kmer_size )

def create_figures( datalist, dir, k ):
    figures = list()
    os.mkdir( dir )
    count = 1
    to_log = lambda x: np.log10( x ) + 1 if x > 0 else 0
    for tup in datalist:
        label, counts = tup
        x = range( 0, len( counts ) )

        fig = plt.figure( count, figsize = ( 5, 4 ) )
        ax1 = fig.add_subplot( 111 )
        ax1.set_title( label + ' (k = %d)' % k, fontsize = 20 )
        ax1.tick_params( which = 'major', labelsize = 15 )
        ax1.set_xlim( 0, len( counts ) )
        ax1.set_ylim( 0, max( [ to_log( x ) for x in counts ] ) + 1 )
        if min( [ to_log( x ) for x in counts ] ) == 1:
            print( label )
        ax1.set_xlabel( 'Position in Sequence', fontsize = 15 ) 
        ax1.set_ylabel( 'log10(Coverage Depth) + 1', fontsize = 15 )
        ax1.plot( x, [ to_log( x ) for x in counts ] )
        count += 1

        plt.savefig( dir + '/' + label + '.pdf', bbox_inches = 'tight' )
        plt.close( fig )
        
def get_positional_counts( kmers, kmer_dict ):
    size = max( [ item.get_start() for item in kmers ] )
    counts = [ 0 ] * size

    for k in kmers:
        try:
            kseq = k._seq

            kmer_count = kmer_dict[ kseq ]
            counts[ k.get_start() - 1 ] = kmer_count
        except KeyError:
            pass
    return counts

def kmer_dict_wcounts( sequences, k ):
    seqs = [ item.get_seq() for item in sequences ]
    kmer_dict = {}

    for seq in seqs:
        kmers = oligo.subset_lists_iter( seq, k, 1 )
        for km in kmers:
            if km not in kmer_dict:
                kmer_dict[ km ] = 0
            kmer_dict[ km ] += 1
    return kmer_dict

def parse_fasta( fname ):
    names, sequences = oligo.read_fasta_lists( fname )
    sequences = [ Sequence( n, s.replace( '-', '' ) ) for n, s in zip( names, sequences ) ] 

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
        end += 1
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
        return self._seq == other._seq and \
            self._start == other._start and \
            self._end == other._end

    def __ne__( self, other ):
        return not( self.__eq__( other ) )

    def __lt__( self, other ):
        return self._start < other._start

if __name__ == '__main__':
    main()
