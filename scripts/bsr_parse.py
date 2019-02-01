#!/usr/bin/env python3
from Bio.Blast import NCBIXML # For parsing blast results
import argparse               # For command-line arguments 

def main():

    arg_parser = argparse.ArgumentParser( description = "Given blast output, determine and output good matches." )

    arg_parser.add_argument( '-b', '--blast_output', help = "File output by BLAST to parse." )


    args = arg_parser.parse_args()

    records = parse_blast( args.blast_output )

    self_scores = find_self_scores( records )
    
def find_self_scores( blast_records ):
    self_hits = set()
    self_hit_names = set()
    num = 0
    for record in blast_records._records:
        for hit_recs in record:
            for hit in hit_recs:
                num += 1
                if self_hit( hit ):
                    add_self_hit( self_hits, hit )
                    self_hit_names.add( hit.query_name )

    return self_hits

def add_self_hit( hit_set, single_hit ):
    new_score = SelfHitScore( name      = single_hit.query_name,
                              hit_score = single_hit.hsp_score
                            )
    hit_set.add( new_score )


def self_hit( hit ):
    split_subject = set( hit.subject_name.split( '|' ) )
    split_query   = set( hit.query_name.split( '|' ) )

    return len( split_subject & split_query )

def parse_blast( blast_file ):

    record_creator = BlastRecordCreator()
    record_creator.create_record( blast_file )

    record = record_creator.as_list()[ 0 ]
    parser = BlastRecordParser()

    return parser.parse( record )
       
class HitScore:
    def __init__( self, name = "",
                  hit_score = 0 
                  ):
        self._name      = name
        self._hit_score = hit_score

    def __eq__( self, other ):
        return self._name == other._name and \
               self._hit_score == other._hit_score
    def __hash__( self ):
        return hash( self._name )

class SelfHitScore( HitScore ):
    def __init__( self, name = "",
                  hit_score = 0
                ):
        super().__init__( name = name, hit_score = hit_score )



class BlastRecordCreator:
    def __init__( self ):
        self._factory = BlastRecordFactory()
        self._blast_data = list()

    def create_record( self, record_name ):
        rec = self._factory.create_record( record_name )
        self._blast_data.append( rec )

    def set_factory( self, new_fac ):
        self._factory = new_factory
    def as_list( self ):
        return self._blast_data
        
class BlastRecordFactory:
    def create_record( self, record_name ):
        return BlastRecord( record_name )

class BlastRecordParser:
    def __init__( self ):
        self._num_hits       = 10
        self._num_hsps       = 1
        self._identity_score = 0.5

    class BlastRecordHit:
        def __init__( self,
                      query_name       = None,
                      query_length     = None,
                      subject_name     = None,
                      subject_length   = None,
                      alignment_length = None,
                      query_start      = None,
                      query_end        = None,
                      subject_start    = None,
                      subject_end      = None,
                      hsp_score        = None,
                      hsp_expect       = None,
                      hsp_identities   = None,
                      percent_match    = None,
                      number_of_gaps   = None
                     ):
            self.query_name       = query_name
            self.query_length     = query_length
            self.subject_name     = subject_name
            self.subject_length   = subject_length
            self.alignment_length = alignment_length
            self.query_start      = query_start
            self.query_end        = query_end
            self.subject_start    = subject_start
            self.subject_end      = subject_end
            self.hsp_score        = hsp_score
            self.hsp_expect       = hsp_expect
            self.hsp_identities   = hsp_identities
            self.percent_match    = percent_match
            self.number_of_gaps   = number_of_gaps

        def __str__( self ):
            
            return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % \
                   ( self.query_name, self.query_length, self.subject_name,
                     self.subject_length, self.alignment_length, self.query_start,
                     self.query_end, self.subject_start, self.subject_end,
                     self.hsp_score, self.hsp_expect, self.hsp_identities, self.percent_match,
                     self.number_of_gaps
                   )

    def parse( self, blast_record ):
        records = list()

        with open( blast_record.get_filename() ) as open_file:
            try:
                
                blast_records = NCBIXML.parse( open_file )

                for record in blast_records:
                    parsed_record = self._parse_record( record )
                    records.append( parsed_record )
            except ValueError as error:
                print( error, blast_record.get_filename() )


        blast_record._records = records

        return blast_record


    def _parse_record( self, record ):
        num_hits = len( record.alignments )
        hits = list()

        if len( record.alignments ) > 0:
            if len( record.alignments ) > self._num_hits:
                num_hits = self._num_hits

            for current_hit in range( num_hits ):
                hits.append( self._parse_hit( record, current_hit ) )

        return hits

    def _parse_hit( self, record, hit ):
        alignment = record.alignments[ hit ]
        num_hsps = min( self._num_hsps, len( alignment.hsps ) )
        out_list = list()

        for this_hsp in range( num_hsps ):
            first_hsp = alignment.hsps[ this_hsp ]

            hit = BlastRecordParser.BlastRecordHit( 
                     query_name       = record.query,
                     query_length     = record.query_length,
                     subject_name     = alignment.title,
                     subject_length   = alignment.length,
                     alignment_length = first_hsp.align_length,
                     query_start      = first_hsp.query_start,
                     query_end        = first_hsp.query_end,
                     subject_start    = first_hsp.sbjct_start,
                     subject_end      = first_hsp.sbjct_end,
                     hsp_score        = first_hsp.score,
                     hsp_expect       = first_hsp.expect,
                     hsp_identities   = first_hsp.identities,
                     percent_match    = float( first_hsp.identities - first_hsp.gaps )/ int( first_hsp.align_length ),
                     number_of_gaps   = first_hsp.gaps
                    )
            out_list.append( hit )
        return out_list

         

    def set_num_hits( self, new_hits ):
        self._num_hits = new_hits

    def set_num_hsps( self, new_hsps ):
        self._num_hsps = new_hsps

    def set_id_score( self, new_id ):
        self._identity_score = new_i


class BlastRecord:
    parser = BlastRecordParser()
    def __init__( self, outfile_name ):
        self._file        = outfile_name
        self._records = None

    def get_filename( self ):
        return self._file

    def get_id( self ):
        return self._file.split( '/' )[ 1 ] \
                   .split( '.' )[ 0 ]

    def set_filename( self, newfile ):
        self._file = newfile

    def get_parser( self ):
        return BlastRecord.parser

    def records_as_string( self ):
        out_string = ""

        for outer_item in self._records:
            for inner_item in outer_item:
                for record in inner_item:
                    out_string += "%s\n" % ( str( record ) )
        return out_string
 

if __name__ == '__main__':
    main()

