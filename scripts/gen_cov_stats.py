#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse
import os
import sys 


def main():
    arg_parser = argparse.ArgumentParser( description = "Parse coverage stats for each cluster processed by a library design." )

    arg_parser.add_argument( '-c', '--cluster_dir', help = "Directory containing original, un-aligned clusters. Each cluster in this directory should have "
                                                           "a corresponding cluster in the design dictionary"
                           )
    arg_parser.add_argument( '-d', '--design_dir', help  = "Directory containing designed clusters, each cluster containing oligos." )
    arg_parser.add_argument( '-o', '--output', help = 'Tab-delimited file to write output to', default = 'output.txt' )

    args = arg_parser.parse_args()

    design_clusters = parse_fasta_dir( args.design_dir )
    orig_clusters   = parse_fasta_dir( args.cluster_dir )

    kmer_sizes = [ 7, 8, 9, 10, 30 ]

    design_stats = get_fasta_stats( design_clusters, kmer_sizes )
    orig_stats   = get_fasta_stats( orig_clusters, kmer_sizes )

    with open( args.output, 'w' ) as open_file:
        header = 'Run Name\tFile Name\tNum 30mers in design\tPerc 7mers in design\tPerc 8mers in design\tPerc 9mers in design\tPerc 10mers in design\tPerc 30mers in design\n'
        open_file.write( header )

        run_name = args.design_dir.strip( '/' )
        for index, stat in enumerate( design_stats ):
            filename = stat.filename
            print( filename )
            print( orig_stats[ index ].filename )
            write_str  = str( stat.num_kmers[ 30 ] ) + '\t'
            write_str += '\t'.join( [ str( item ) for item in stat.divide( orig_stats[ index ] ) ] )
            open_file.write( run_name + '\t' + filename + '\t' + write_str + '\n' )
        



class Sequence:
    def __init__( self, name = "", sequence = "" ):
        self.name     = name
        self.sequence = sequence
    def get_kmers_with_step( self, k, step_size ):
        return oligo.subset_lists_iter( self.sequence, k, step_size )

    def __hash__( self ):
        return hash( self.sequence )
    def __eq__( self, other ):
        return self.sequence == other.sequence
    def __ne__( self, other ):
        return not self.__eq__( other )
    def __str__( self ):
        return '>%s\n%s\n' % ( self.name, self.sequence )
    def __len__( self ):
        return len( self.sequence )
    
class FastaParser:
    instance = None
    def __init__( self ):
        pass

    def get_instance():
        if FastaParser.instance is None:
            FastaParser.instance = FastaParser()
        return FastaParser.instance
        
    def parse( self, filename ):
        out_fasta = Fasta( filename = filename )

        names, sequences = oligo.read_fasta_lists( filename )
        for name, seq in zip( names, sequences ):
            out_fasta.add_seq( Sequence( name = name, sequence = seq ) )
        return out_fasta

class Fasta:
    def __init__( self, filename = "", sequences = None ):
        self.filename = filename.split( '/' )[ -1 ]
        self.sequences = list() if sequences is None else sequences

    def add_seq( self, new_sequence ):
        self.sequences.append( new_sequence )

    def get_sequences( self ):
        return self.sequences

class FastaStats:
    def __init__( self, filename ):
        self.filename  = filename
        self.num_kmers = dict()

    def __hash__( self ):
        return hash( filename )

    def __eq__( self, other ):
        return self.filename == other.filename 
    def divide( self, other ):
        out_vals = list()

        for key, value in self.num_kmers.items():
            out_vals.append( value / other.num_kmers[ key ] )
        return out_vals
            

    def __str__( self ):
        out_str = '%s\t' % self.filename

        for k, count in self.num_kmers.items():
            out_str += '%d\t%d' % ( k, count )
        return out_str + '\n'

def parse_fasta_dir( dirname ):
    files      = sorted( os.listdir( dirname ) )
    parser     = FastaParser.get_instance()
    out_fastas = list()

    for file in files:
        out_fastas.append( parser.parse( dirname + '/' + file ) )
    return out_fastas
        
        
def get_fasta_stats( clusters, kmer_sizes ):
    fasta_stats = list()
    for fasta in clusters:
        fasta_stats.append( get_fasta_stat( fasta, kmer_sizes ) )
    return fasta_stats

def get_fasta_stat( fasta, kmer_sizes ):
    stat = FastaStats( fasta.filename )
    seqs = fasta.get_sequences()

    for k in kmer_sizes:
        k_set = set()

        for seq in seqs:
            k_set |= seq.get_kmers_with_step( k, 1 )
        stat.num_kmers[ k ] = len( k_set )
    return stat
            
    

if __name__ == '__main__':
    main()
