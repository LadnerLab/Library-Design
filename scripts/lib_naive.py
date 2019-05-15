#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse


def main():
    arg_parser = argparse.ArgumentParser( description = "Simple creation of a library, both gap spanning and non-gap spanning algorithms" )

    arg_parser.add_argument( '-q', '--query', help = "Fasta query file." )
    arg_parser.add_argument( '-o', '--output', help = "Fasta file to output" )
    arg_parser.add_argument( '-g', '--gap_span', help = "Fasta query file.", default = False, action = "store_true" )
    arg_parser.add_argument( '-w', '--window_size', help = "Window Size to use for grabbing oligos.", default = 19, type = int )
    arg_parser.add_argument( '-s', '--step_size', help = "Number of amino acids to step after each window.", default = 10, type = int )

    args = arg_parser.parse_args()

    names, sequences = oligo.read_fasta_lists( args.query )
    seqs = list()

    for name, sequence in zip( names, sequences ):
        seqs.append( Sequence( name = name, sequence = sequence ) )

    print( "Number of input sequences: ", len( seqs ) )

    if args.gap_span:
        designer = GapSpanningLibraryDesigner( window_size = args.window_size, step_size = args.step_size )
    else:
        designer = LibraryDesigner( window_size = args.window_size, step_size = args.step_size )

    library = designer.design( seqs )

    print( "Number of output Kmers: ", len( library ) )

    write_sequences( args.output, library )

class LibraryDesigner():
    def __init__( self, window_size = 0, step_size = 0 ):
        self.window_size = window_size
        self.step_size   = step_size

    def _get_oligos( self, sequence ):
        xmers = set()

        start = 0
        end = self.window_size

        seq = sequence.sequence

        if len( seq ) < self.window_size and 'X' not in seq:
            xmers.add( seq )
        while end < len( seq ):
            xmer = seq[ start:end ]

            new_name = sequence.name + "_%d_%d" % ( start, end )
            if not 'X' in xmer and '-' not in xmer:
                xmers.add( Sequence( name = new_name,
                                     sequence = xmer
                                   )
                         )
            start += self.step_size
            end   = start + self.window_size
        return xmers

    def design( self, sequences ):
        all_oligos = set()

        for seq in sequences:
            oligos = self._get_oligos( seq )
            all_oligos |= oligos
        return all_oligos

    def _valid_oligo( oligo ):
        return not 'X' in oligo
        
class GapSpanningLibraryDesigner( LibraryDesigner ):
    def __init__( self, window_size = 0, step_size = 0 ):
        super().__init__( window_size = window_size, step_size = step_size )

    def design( self, sequences ):
        all_oligos = set()

        for seq in sequences:
            oligos = self._get_oligos( seq )
            all_oligos |= oligos
        return all_oligos

    def _get_oligos( self, seq ):
        start = 0
        sequence = seq.sequence
        oligos = set()

        while start + self.window_size <= len( seq ):
            current   = start
            probe     = start
            cur_oligo = ""

            # current - start = size of oligo
            while current - start < self.window_size and probe < len( seq ):
                if sequence[ probe ] != '-':
                    cur_oligo += sequence[ probe ]
                    current   += 1
                probe += 1

            if len( cur_oligo ) == self.window_size and 'X' not in cur_oligo:
                new_name = seq.name + "_%d_%d" % ( start, probe )
                oligos.add( Sequence( name = new_name, sequence = cur_oligo ) )

            start += self.step_size
        return oligos
            


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
    

def write_sequences( filename, library ):
    with open( filename, 'w' ) as open_file:
        for seq in library:
            open_file.write( str( seq ) )


if __name__ == '__main__':
    main()
