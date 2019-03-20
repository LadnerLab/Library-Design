#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse

def main():
    arg_parser = argparse.ArgumentParser( description = "Verify that the positions of oligos in a "
                                                        "designed library actually appear at the "
                                                        "positions at which they are reported."
                                        )
    arg_parser.add_argument( '-o', '--original_seqs', help = "Sequences that were used to "
                                                             "create the oligo library"
                           )
    arg_parser.add_argument( '-l', '--library', help = "Designed oligonucleotide library "
                                                       "containing either encoded nt sequences, "
                                                       "or amino acid sequences. Note that if nt "
                                                       "sequences are included, a map translating "
                                                       "amino acids to codons and a map translating "
                                                       "encoded oligo names to un-encoded oligo names "
                                                       "must also be provided."
                           )
    arg_parser.add_argument( '-a', '--acid_map', help = "Map translating amino acids to "
                                                        "nucleotide codons. Must be in the form:  "
                                                        "{letter},{nucleotides,3},{weighting},{index}. "
                                                        "The weightings do not need to sum to 1."
                                                        "Codon indices must range from 0 to 63."
                           )
       

    args = arg_parser.parse_args()
    fasta_parse = FastaParser()

    oligo_seqs    = fasta_parse.parse( args.library )
    original_seqs = fasta_parse.parse( args.original_seqs )
    loc_names = SequenceWithLocation.add_locs_to_seq_list( oligo_seqs )

class Parser:
    def parse( self, filename ):
        with open( filename, 'r' ) as open_file:
            return open_file.readlines()

class FastaParser( Parser ):
    def parse( self, filename ):
        names, sequences = oligo.read_fasta_lists( filename )
        out_seqs = list()

        for index, name in enumerate( names ):
            out_seqs.append( Sequence( name = name,
                                       sequence = sequences[ index ]
                                     )
                           )
        return out_seqs
            
class Sequence:
    def __init__( self, name = "", sequence = "" ):
        self.name     = name
        self.sequence = sequence

    def __str__( self ):
        return ">%s\n%s" % ( self.name, self.sequence )

class SequenceWithLocation( Sequence ):
    def __init__( self, name = "", sequence = "",
                  location_start = 0, location_end = 0 ):
        self.name     = name
        self.sequence = sequence

        self.location_start = location_start

        if location_end:
            self.location_end = location_end
        else:
            self.location_end = len( sequence )
    def add_locs_to_seq_list( seq_list ):
        out_seqs = list()
        for sequence in seq_list:
            name       = sequence.name
            split_name = name.split( '_' )
            sequence = sequence.sequence
            loc_start = split_name[ -2 ]
            loc_end   = split_name[ -1 ]

            try:
                out_seqs.append( SequenceWithLocation( name = name, sequence = sequence,
                                                       location_start = int( loc_start ),
                                                       location_end   = int( loc_end )
                                                     )
                               )
            except ValueError:
                print( "WARNING: %s does not contain the location of an oligo, and will therefore not "
                       "be considered." % ( name )
                     )
        return out_seqs 

    def __str__( self ):
        return ">%s_%d_%d\n%s" % ( self.name, self.location_start,
                                   self.location_end, self.sequence
                                 )



                
if __name__ == '__main__':
    main()
