#!/usr/bin/env python3
from Bio.Blast import NCBIXML # For parsing blast results
import argparse               # For command-line arguments 
import re                     # for finding nc-tags and tax ids

def main():

    arg_parser = argparse.ArgumentParser( description = "Given blast output, determine and output good matches." )

    arg_parser.add_argument( '-s', '--self_blast', help = "Self BLAST output to parse" )
    arg_parser.add_argument( '-r', '--ref_blast',  help  = "Reference BLAST output to parse" )
    arg_parser.add_argument( '-g', '--good_hit',   help   = "Floating point ratio determining what "
                                                            "a good hit is.",
                             type = float, default = 0.7
                           )
    arg_parser.add_argument( '-n', '--num_hits', help = "Number of hits to consider for each reference blast",
                             type = int, default = 10
                           )
    arg_parser.add_argument( '--rep_to_tax', help = 'Name of file containing NC_tag to tax id mapping, '
                                                    'this file should contain a mapping for each nc_tag found in '
                                                    'the self_blast sequences'
                           )


    args = arg_parser.parse_args()

    IDENTITY = args.good_hit

    nc_taxid     = parse_nc_taxid( args.rep_to_tax )
    self_records = parse_blast( args.self_blast )
    ref_records  = parse_blast( args.ref_blast, num_hits = args.num_hits )

    self_scores     = find_self_scores( self_records )
    non_self_scores = find_good_hits( ref_records, IDENTITY ) 

    non_self_scores_d = scores_to_dict( non_self_scores )

    bsr_hits = get_good_bsr_scores( self_scores, non_self_scores_d, IDENTITY )

    for hit in bsr_hits:
        print( str( hit ) )

    label_bsr_hits_with_taxids( nc_taxid, bsr_hits )


def label_bsr_hits_with_taxids( nc_taxid, bsr_hits ):
    query_pattern = r'OXX=[0-9]*,[0-9]*,[0-9]*,[0-9]*'

    suffix = '[0-9]+\.[0-9]'
    ref_patterns = ( r'NC_%s' % suffix,
                     r'KJ%s'  % suffix,
                     r'LC%s'  % suffix,
                     r'MG%s'  % suffix,
                     r'AB%s'  % suffix,
                     r'KY%s'  % suffix,
                     r'LT%s'  % suffix
                   )

    for hit in bsr_hits:
        query = hit._query
        ref   = hit._ref

        tax_ids = re.search( query_pattern, query ).group()

        query_id = get_id_from_string( tax_ids )

def get_id_from_string( string ):
    ids = string.split( 'OXX=' )[ 1 ]

    id_list = ids.split( ',' )

    # return species id, if available
    if id_list[ 1 ]:
        return id_list[ 1 ]

    # return the first id
    return id_list[ 0 ]

def get_good_bsr_scores( self_score_set, non_self_dict, good_hit_thresh ):
    out_list = list()
    for score in self_score_set:
        if score in non_self_dict:
            for current in non_self_dict[ score ]:
                bsr = calc_bsr( current, score )
                if bsr >= good_hit_thresh:
                    out_list.append( BSRScore( query = current._name,
                                               ref = score._name,
                                               bsr = bsr
                                             )
                                    )
    return out_list

class BSRScore:
    def __init__( self, query = None, ref = None,
                  bsr = None
                ):
        self._query = query
        self._ref   = ref
        self._bsr   = bsr

        self._query_id = None
        self._ref_id   = None

    def set_query_id( self, id ):
        self._query_id = id

    def set_ref_id( self, id ):
        self._ref_id = id 

    def __str__( self ):
        return '%s\t%s\t%f' % ( self._query, self._ref, self._bsr )

def parse_nc_taxid( filename ):
    out_dict = {}
    with open( filename, 'r' ) as open_file:
        for lineno, line in enumerate( open_file ):
            if lineno:
                nc_tag, taxid = line.strip().split( '\t' )
                out_dict[ nc_tag ] = taxid
    return out_dict

def calc_bsr( query, ref ):
    return query._hit_score / ref._hit_score 

def scores_to_dict( list_of_scores ):
    out_dict = {}
    scores = list_of_scores

    for item in scores:
        if item not in out_dict:
            out_dict[ item ] = list()
        out_dict[ item ].append( item )

    return out_dict


def find_good_hits( ref_recs, identity_score ):
    hit_scores = list()

    for record in ref_recs._records:
        for hits in record:
            for hit in hits:
                hit_scores.append( add_good_hit( hit ) ) # Add all hits, parse good ones out later
    return hit_scores

def good_hit( hit, identity_score ):
    return hit.percent_match >= identity_score

def add_good_hit( hit ):
    return HitScore( name = hit.query_name,
                     other_name = hit.subject_name,
                     hit_score = hit.hsp_score
                   )

def find_self_scores( blast_records ):
    self_hits = set()
    for record in blast_records._records:
        for hit_recs in record:
            for hit in hit_recs:
                if self_hit( hit ):
                    add_self_hit( self_hits, hit )

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

def parse_blast( blast_file, num_hits = None ):
    record_creator = BlastRecordCreator()
    record_creator.create_record( blast_file )

    record = record_creator.as_list()[ 0 ]
    parser = BlastRecordParser()

    if num_hits:
        parser.set_num_hits( num_hits )
    return parser.parse( record )
       
class HitScore:
    def __init__( self, name = "",
                  other_name = "",
                  hit_score = 0 
                  ):
        self._name       = name
        self._other_name = other_name
        self._hit_score  = hit_score

    def __eq__( self, other ):
        return hash( self ) == hash( other )
    def __hash__( self ):
        if '|' in self._other_name:
            return hash( self._other_name.split( '|' )[ -1 ] )
        return hash( self._other_name )
    
    def __str__( self ):
        return '%s\t%s\t%d' % ( self._name, self._other_name,
                                self._hit_score
                              )
    def __contains__( self, key ):
        return key in self._other_name

class SelfHitScore( HitScore ):
    def __init__( self, name = "",
                  hit_score = 0
                ):
        super().__init__( name = name, hit_score = hit_score )

    def __hash__( self ):
        return hash( self._name.split( '|' )[ -1 ] )

    def __str__( self ):
        return '%s\t%d' % ( self._name, self._hit_score )
    def __contains__( self, key ):
        return key in self._name

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
        self._identity_score = new_id


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

