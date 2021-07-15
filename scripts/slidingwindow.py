#!/usr/bin/env python3

import fastatools as ft  # Available here: https://github.com/jtladner/Modules
import argparse


def main():
    arg_parser = argparse.ArgumentParser( description = "Peptide design using a sliding window approach.", formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    arg_parser.add_argument( '-w', '--window_size', help = "Length of desired peptides.", default = 30, type = int )
    arg_parser.add_argument( '-s', '--step_size', help = "Number of amino acids to move between each window.", default = 1, type = int )
    arg_parser.add_argument( '-g', '--gap_span', help = "Use this flag if you want to use the gap-spanning approach for peptide design.", default = False, action = "store_true" )
    arg_parser.add_argument( '-q', '--quiet', help = "Use this flag if you do not want any info printed to screen during run time.", default = False, action = "store_true" )

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-t', '--targets', help = "Fasta file contianing target protein sequences.", required=True )
    reqArgs.add_argument( '-o', '--output', help = "Name for output file containing designed peptides (fasta format).", required=True )

    args = arg_parser.parse_args()

    names, sequences = ft.read_fasta_lists( args.targets )
    seqs = list()

    for name, sequence in zip( names, sequences ):
        seqs.append( Sequence( name = name, sequence = sequence ) )

    if not args.quiet:
        print( "Number of input sequences: ", len( seqs ) )

    if args.gap_span:
        designer = GapSpanningLibraryDesigner( window_size = args.window_size, step_size = args.step_size )
    else:
        designer = LibraryDesigner( window_size = args.window_size, step_size = args.step_size )

    library = designer.design( seqs )

    if not args.quiet:
        print( "Number of output Kmers: ", len( library ) )

    outD = {e.name:e.sequence for e in library}
    namesSorted = sorted(list(outD.keys()))
    ft.write_fasta(namesSorted, [outD[n] for n in namesSorted], args.output)

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

            new_name = sequence.name + "_%03d_%03d" % ( start, end )
            if not 'X' in xmer and '-' not in xmer:
                xmers.add( Sequence( name = new_name, sequence = xmer))
                
                # If the remaining sequence is less than the step size
                # Add a peptide covering the end of the seq
                if  0 < len( seq[end:] ) < (self.step_size):
                    xmer = seq[-self.window_size:]
                    new_name = sequence.name + "_%03d_%03d" % ( len(seq)-self.window_size, len(seq) )
                    xmers.add( Sequence( name = new_name, sequence = xmer))
                
            start += self.step_size
            end   = start + self.window_size
        return xmers

    def design( self, sequences ):
        all_oligos = set()

        for seq in sequences:
            oligos = self._get_oligos( seq )
            all_oligos |= oligos
        return all_oligos

        
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

        while start + self.window_size <= len( seq.sequence ):
            current   = start
            probe     = start
            cur_oligo = ""

            # current - start = size of oligo
            while current - start < self.window_size and probe < len( seq.sequence ):
                if sequence[ probe ] != '-':
                    cur_oligo += sequence[ probe ]
                    current   += 1
                probe += 1

            if len( cur_oligo ) == ( self.window_size ) and 'X' not in cur_oligo:

                # If the remaining sequence is less than the step size
                if 0 < len( sequence[ probe: ].replace( '-', '' ) ) < (self.step_size):

                    # Go ahead and add the current peptide
                    new_name = seq.name + "_%03d_%03d" % ( start, probe )
                    oligos.add( Sequence( name = new_name, sequence  = cur_oligo ) )
                    
                    # Adjust variables to add a final peptide with a larger than typical overlap
                    cur_oligo += sequence[ probe: ].replace( '-', '' )
                    start += len( cur_oligo ) - self.window_size  
                    cur_oligo = cur_oligo[ len( cur_oligo ) - self.window_size: ]
                    probe = start + len(cur_oligo)

                new_name = seq.name + "_%03d_%03d" % ( start, probe )
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
    


if __name__ == '__main__':
    main()
