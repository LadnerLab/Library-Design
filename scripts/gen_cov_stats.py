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

    args = arg_parser.parse_args()

    design_clusters = parse_fasta_dir( args.design_dir )

class Sequence:
    def __init__( self, name = "", sequence = "" ):
        self.name     = name
        self.sequence = sequence
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
        self.filename = ""
        self.sequences = list() if sequences is None else sequences

    def add_seq( self, new_sequence ):
        self.sequences.append( new_sequence )

def parse_fasta_dir( dirname ):
    files      = os.listdir( dirname )
    parser     = FastaParser.get_instance()
    out_fastas = list()

    for file in files:
        out_fastas.append( parser.parse( dirname + '/' + file ) )
    return out_fastas
        
        

if __name__ == '__main__':
    main()
