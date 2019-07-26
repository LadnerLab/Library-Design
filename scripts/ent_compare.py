#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "" )

    arg_parser.add_argument( '-k', '--kmer_sizes', help = "The kmer size to use for evaluations." )
    arg_parser.add_argument( '-d', '--designed', help = "Fasta file containing designed sequences." )
    arg_parser.add_argument( '-c', '--check', help = "Fasta file containing sequence to check coverage." )
    arg_parser.add_argument( '--create_table', action = 'store_true' )
    arg_parser.add_argument( '-o', '--output', default = "out.txt" )

    args = arg_parser.parse_args()

    # if a table is created choose to do that
    if args.create_table:
        output = create_table( args )
        write_kmer_counts( args.output, output, args.kmer_sizes )


def create_table( args ):
    k_sizes = args.kmer_sizes.split( ',' )
    to_int = lambda x: [ int( i ) for i in x ] 
    k_sizes = to_int( k_sizes )


    check_file  = args.check
    design_file = args.designed

    designed_sequences = parse_fasta( args.designed )
    check_sequences    = parse_fasta( args.check )

    kmer_counts = dict()

    for k in k_sizes:
        designed_kmers = designed_sequences.apply( get_kmers, k )

        designed_all = set()
        check_all    = set()
        for kmers in designed_kmers:
            designed_all |= kmers

        for seq in check_sequences:
            kmers = get_kmers( seq, k )
            covered = kmers & designed_all
            uncovered = len( kmers - covered )

            if seq.name not in kmer_counts:
                kmer_counts[ seq.name ] = dict()
            kmer_counts[ seq.name ][ k ] = ( len( covered ), len( kmers ) )
    return kmer_counts
    
def write_kmer_counts( fname, output_dict, k_sizes ):
    with open( fname, 'w' ) as open_file:
        header = 'Sequence Name ' + '\t'
        k_sizes = sorted( [ int( item ) for item in k_sizes.split( ',' ) ] )
        k_str   = '\t'.join( [ 'Percent ' + str( item ) + 'mers covered by design' for item in k_sizes ] )
        header += k_str + '\n'
        open_file.write( header )

        to_str = lambda x: str( 100 * ( x[ 0 ] / x[ 1 ] ) )

        for seq_name, kmer_stats in output_dict.items():

            open_file.write( seq_name + '\t' +
                             '\t'.join( [ to_str( kmer_stats[ item ] ) for item in sorted( list( kmer_stats.keys() ) ) ] )
                             + '\n'
                           )



def get_kmers( seq, k ):
    s = seq.seq
    return oligo.subset_lists_iter( s, k, 1 )
    
def parse_fasta( fname ):
    n, s = oligo.read_fasta_lists( fname )
    seqs = [ Sequence( name = na, seq = se.replace( '-', '' ) ) for na, se in zip( n, s ) ] 
    return SequenceCollection( seqs )

class SequenceCollection:
    def __init__( self, sequences ):
        self.seqs = sequences

    def apply( self, func, *args, **kwargs ):
        return [ func( i, *args, **kwargs ) for i in self.seqs ] 

    def __iter__( self ):
        return iter( self.seqs )
    def __next__( self ):
        return next( self.seqs )
    
class Sequence:
    def __init__( self, name = "", seq = "" ):
        self.name = name
        self.seq  = seq

if __name__ == '__main__':
    main()
