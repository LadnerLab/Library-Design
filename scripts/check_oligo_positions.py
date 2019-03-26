#!/usr/bin/env python3
import protein_oligo_library as oligo
import argparse
import re

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
    arg_parser.add_argument( '-e', '--encoding_map', help = "Map translating encoded names into the "
                                                            "names of sequences which have been encoded."
                           )
       

    args = arg_parser.parse_args()
    fasta_parse = FastaParser()

    if ( args.acid_map and args.encoding_map ):
        aa_map     = parse_acid_map( args.acid_map )
        parsed_name_map = parse_name_map( args.encoding_map )
        oligo_seqs = encoding_to_aa( csv_to_seqs( args.library ), aa_map )
        oligo_seqs = fix_seq_names( oligo_seqs, parsed_name_map )
    else:
        oligo_seqs    = fasta_parse.parse( args.library )
        
    original_seqs = fasta_parse.parse( args.original_seqs )
    loc_names     = SequenceWithLocation.add_locs_to_seq_list( oligo_seqs )
    orig_seqs_dict = seqs_to_dict( original_seqs )

    bad_seqs = validate_oligo_locations( orig_seqs_dict, loc_names )

    if len( bad_seqs ) > 0:
        report_bad_seqs( bad_seqs )
    else:
        print( "No errors detected." )

def fix_seq_names( seqs, name_map ):
    out_seqs = list()
    for sequence in seqs:
        orig_name = sequence.name
        new_name  = name_map[ orig_name ]
        out_seqs.append( Sequence( name = new_name,
                                   sequence = sequence.sequence
                                 )
                       )
    return out_seqs
        
def parse_name_map( filename ):
    out_dict = {}
    with open( filename, 'r' ) as open_file:
        for line in open_file:
            split_line = line.split( '\t' )
            encod_name = split_line[ 1 ]
            orig_name    = split_line[ 0 ]

            out_dict[ encod_name ] = orig_name
    return out_dict

def encoding_to_aa( seqs, aa_map ):
    out_seqs = list()
    for seq in seqs:
        aa_seq = nt_to_aa( seq.sequence, aa_map )
        new_sequence = Sequence( name = seq.name, sequence = aa_seq )
    return out_seqs

def nt_to_aa( sequence, aa_map ):
    out_str = ""
    aa_per_codon = 3

    for index in range( 0, len( sequence ) // aa_per_codon, aa_per_codon ):
        codon = sequence[ index : index + aa_per_codon ]
        out_str += aa_map[ codon ]
    return out_str
    
def csv_to_seqs( filename ):
    with open( filename, 'r' ) as open_file:
        out_seqs = list()
        for lineno, line in enumerate( open_file ):
            if lineno: # skip header
                split_line = line.split( ',' )
                name = split_line[ 0 ].strip()
                nt_encoding = split_line[ 2 ].strip()
                name = name.replace( '~', '_' )
                name = ''.join( name.split( '_' )[ :-1 ] )

                out_seqs.append( Sequence( name = name, sequence = nt_encoding ) )
        return out_seqs
    

def parse_acid_map( filename ):
    out_map = {}
    with open( filename, 'r' ) as open_file:
        for line in open_file:
            split_line = line.split( ',' )
            codon = split_line[ 1 ]
            aa    = split_line[ 0 ]
            out_map[ codon ] = aa 
    return out_map

def report_bad_seqs( seq_list ):
    print( "Original Name\tOligo Name\tOriginal Seq\tOligo" )
    for item in seq:
        print( "%s\t%s\t%s\t%s\t" % ( item.original_seq.name,
                                      item.oligo_seq.name,
                                      item.origina._seq.sequence,
                                      item.oligo_seq.sequence
                                    )
             )
def seqs_to_dict( sequence_list ):
    out_dict = {}
    for seq in sequence_list:
        out_dict[ seq.name ] = seq
    return out_dict

def validate_oligo_locations( originals, with_locs ):
    error_seqs = list()
    for sequence in with_locs:
        loc_start = sequence.location_start
        loc_end   = sequence.location_end

        oligo     = sequence.sequence

        try:
            orig_seq_oligo = originals[ sequence.name ].sequence[ loc_start:loc_end ]
            assert( oligo == orig_seq_oligo )
        except AssertionError:
            error_seqs.append( ErrorSequence( original_seq = Sequence( name = sequence.name,
                                                                       sequence = orig_seq_oligo
                                                                     ),
                                              oligo_seq    = Sequence( name = sequence.name_str(),
                                                                       sequence = oligo
                                                                     )
                                            )
                             )
        except KeyError:
            pass
    return error_seqs

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

class ErrorSequence( Sequence ):
    def __init__( self, original_seq = "", oligo_seq = "" ):
        self.original_seq = ""
        self.oligo_seq    = ""

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
    def name_str( self ):
               out_string = "%s_%d_%d" % ( self.name,
                                           self.location_start,
                                           self.location_end
                                         )
               return out_string

    def add_locs_to_seq_list( seq_list ):
        pattern = re.compile ( '_(\d+)_(\d+)' )
        out_seqs = list()
        for sequence in seq_list:
            name       = sequence.name
            split_name = pattern.search( name )
            sequence = sequence.sequence
            loc_start = split_name.group( 1 )
            loc_end   = split_name.group( 2 )
            new_name  = name.split( split_name.group() )[ 0 ]

            try:
                out_seqs.append( SequenceWithLocation( new_name, sequence = sequence,
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
