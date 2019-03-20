#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse

def main():
    arg_parser = argparse.ArgumentParser( description = "Verify that the positions of oligos in a "
                                                        "designed library actually appear at the "
                                                        "positions at which they are reported."
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

    oligo_seqs = fasta_parse.parse( args.library )

class FastaParser:
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

                
if __name__ == '__main__':
    main()
