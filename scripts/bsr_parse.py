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
    arg_parser.add_argument( '--invert', help = "Include flag if you want to invert the good hit criteria. By default, "
                                                "a good hit has a BSR of >= good_hit, inclusion of this flag sets this test to "
                                                " BSR < good_hit", default = False, action = 'store_true'
                           )
    arg_parser.add_argument( '-o', '--output', help = "File to write output to, output format is of the form:"
                                                      "re_name design_name refseq_taxid design_taxid bsr",
                             default = 'bsr_output.tsv' 
                           )
    arg_parser.add_argument( '-t', '--taxonomic_level', help = "Taxonomic level at which to compare hits. "
                                                              "Must be one of either 'family', 'genus', or 'species'.",
                            default = 'species'
                          ) 
     


    args = arg_parser.parse_args()

    IDENTITY = args.good_hit

    tax_rank     = get_tax_rank( args.taxonomic_level )

    nc_taxid     = parse_nc_taxid( args.rep_to_tax, tax_rank )
    self_records = parse_blast( args.self_blast, num_hits = args.num_hits )
    ref_records  = parse_blast( args.ref_blast, num_hits = args.num_hits )

    self_scores     = find_self_scores( self_records )
    non_self_scores = find_good_hits( ref_records, IDENTITY ) 

    non_self_scores_d = scores_to_dict( non_self_scores ) # removes duplicates

    best_hits = get_best_hits( non_self_scores_d )

    bsr_hits = get_good_bsr_scores( self_scores, best_hits, IDENTITY, inverted = args.invert )

    hits_with_ratio = label_bsr_hits_with_taxids( nc_taxid, bsr_hits, tax_rank )

    mismatch_ids = get_hits_mismatch_taxids( hits_with_ratio )

    write_biggest_hits( mismatch_ids, args.output )


def get_tax_rank( level_str ):
    levels = { 'species': 1,
               'genus':   2,
               'family':  3
             }
    level = level_str.lower()

    return levels[ level ]

def get_best_hits( hit_dict ):
    out_dict = {}

    for name, hits in hit_dict.items():
        best_score = 0
        new_hits = list()

        for hit in hits:
            if hit._hit_score > best_score:
                new_hits = list()
                new_hits.append( hit )
                best_score = hit._hit_score
            elif hit._hit_score == best_score:
                new_hits.append( hit )

        out_dict[ name ] = new_hits
                
    return out_dict

def write_biggest_hits( bsr_reports, outfile_name ):
    HEADER = 'Refseq Name\tDesign Name\tRefseq TaxID\t Design TaxID\tBlast Score Ratio\t' + \
             'Percent RefSeq in Alignment\tPercent Design Seq in Alignment\tPercent match'

    with open( outfile_name, 'w' ) as open_file:
        open_file.write( '%s\n' % HEADER )

        for hit in bsr_reports:
            if ( hit._bsr._perc_s_in_align > 0.80 or \
                 hit._bsr._perc_q_in_align > 0.80 ) and \
                 hit._bsr._perc_match > 0.95:

                open_file.write( '%s\n' % str( hit ) )


def get_biggest_hits( bsr_items ):
    out_dict = {}
    for current in bsr_items:
        if current._bsr._query not in out_dict:
            out_dict[ current._bsr._query ] = current._bsr._bsr
        else:
            if out_dict[ current._bsr._query ] < current._bsr._bsr:
                out_dict[ current._bsr._query ] = current._bsr._bsr 
    return out_dict

def get_hits_mismatch_taxids( hits ):
    out_list = list()
    for hit in hits:
        if hit._query_id.strip() != hit._ref_id.strip():
            out_list.append( hit )

    return out_list
        
def label_bsr_hits_with_taxids( nc_taxid, bsr_hits, taxonomic_rank ):
    query_pattern = r'OXX=[0-9]* *,[0-9]* *,[0-9]* *,[0-9]* *'
    out_labelled_hits = list()

    suffix = '[0-9]+\.[0-9]*'
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

        try:
            tax_ids = re.search( query_pattern, query ).group()

            query_id = get_id_from_string( tax_ids, taxonomic_rank )

            ref_id_tag = get_id_tag_from_string( ref, ref_patterns )

            out_labelled_hits.append( BSRScoreWithTaxID( query_id = query_id,
                                                         ref_id   = nc_taxid[ ref_id_tag ],
                                                         bsr = hit
                                                       )
                                 )
        except:
            print( query )
    return out_labelled_hits

class BSRScoreWithTaxID:
    def __init__( self, query_id = None, ref_id = None, bsr = None ):
        self._query_id = query_id
        self._ref_id   = ref_id
        self._bsr      = bsr

    def __hash__( self ):
        return hash( self._bsr ) * hash( self._query_id ) * hash( self._ref_id )
    def __eq__( self, other ):
        return self._query_id == other._query_id and \
               self._ref_id == other._ref_id and \
               self._bsr == other._bsr

    def __str__( self ):
        return '%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f' % (
               self._bsr._ref, self._bsr._query,
               self._ref_id, self._query_id,
               self._bsr._bsr,
               self._bsr._perc_s_in_align * 100,
               self._bsr._perc_q_in_align * 100,
               self._bsr._perc_match      * 100
            )

    def __ne__( self, other ):
        return not self.__eq__( other )

def get_id_tag_from_string( string, patterns ):

    for pattern in patterns:
        matched_str = re.search( pattern, string )
        if re.search( pattern, string ):
            return matched_str.group().split( '.' )[ 0 ]
    return None

def get_id_from_string( string, level ):
    ids = string.split( 'OXX=' )[ 1 ]

    id_list = ids.split( ',' )

    # return species id, if available
    if id_list[ level ]:
        return id_list[ level ]

    # return the first id
    return id_list[ 0 ]

def get_good_bsr_scores( self_score_set, non_self_dict, good_hit_thresh, inverted = False ):
    out_list = list()
    for score in self_score_set:
        for key in non_self_dict.keys():

            if match( score._name, non_self_dict[ key ] ):
                for current in non_self_dict[ key ]:
                    bsr = calc_bsr( current, score )
                    if not inverted and bsr >= good_hit_thresh:
                        out_list.append( BSRScore( query = current._other_name,
                                                   ref = score._name,
                                                   bsr = bsr,
                                                   perc_subject_in_align = current._perc_s_in_align,
                                                   perc_query_in_align   = current._perc_q_in_align,
                                                   percent_match         = current._perc_match
                                                 )
                                        )
                    elif inverted and bsr < good_hit_thresh:
                        out_list.append( BSRScore( query = current._other_name,
                                                   ref = score._name,
                                                   bsr = bsr,
                                                   perc_subject_in_align = current._perc_s_in_align,
                                                   perc_query_in_align   = current._perc_q_in_align,
                                                   percent_match         = current._perc_match
                                                 )
                                        )
                        
    return out_list

def match( name, values ):
    for val in values:
        if name in val._name:
            return True
    return False
    
class BSRScore:
    def __init__( self, query = None, ref = None,
                  bsr = None,
                  perc_query_in_align   = None,
                  perc_subject_in_align = None,
                  percent_match         = None
                ):
        self._query = query
        self._ref   = ref
        self._bsr   = bsr
        self._perc_q_in_align = perc_query_in_align
        self._perc_s_in_align = perc_subject_in_align
        self._perc_match      = percent_match

    def __eq__( self, other ):
        return self._query == other._query and \
               self._ref   == other._ref and \
               self._bsr   == other._bsr

    def __ne__( self, other ):
        return not self.__eq__( other )

    def __hash__( self ):
        return hash( self._query ) * hash( self._ref ) 
               
    def __str__( self ):
        return '%s\t%s\t%f' % ( self._query, self._ref, self._bsr )

def parse_nc_taxid( filename, level ):
    out_dict = {}
    with open( filename, 'r' ) as open_file:
        for lineno, line in enumerate( open_file ):
            if lineno:
                split_line = line.strip().split( '\t' )
                nc_tag = split_line[ 0 ]
                out_dict[ nc_tag ] = ( split_line[ level ] )
    return out_dict

def calc_bsr( query, ref ):
    return query._hit_score / ref._hit_score 

def scores_to_dict( list_of_scores ):
    out_dict = {}
    scores = list_of_scores

    for item in scores:
        if item._other_name not in out_dict:
            out_dict[ item._other_name ] = set()
        out_dict[ item._other_name ].add( item )

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
    align_len = hit.alignment_length
    return HitScore( other_name = hit.query_name,
                     name = hit.subject_name,
                     hit_score = hit.hsp_score,
                     perc_query_in_align   = align_len / hit.query_length,
                     perc_subject_in_align = align_len / hit.subject_length,
                     percent_match         = hit.percent_match
                   )

def find_self_scores( blast_records ):
    self_hits = list()
    for record in blast_records._records:
        for hit_recs in record:
            for hit in hit_recs:
                if self_hit( hit ):
                    add_self_hit( self_hits, hit )

    return self_hits

def add_self_hit( hit_set, single_hit ):
    align_len = single_hit.alignment_length
    new_score = SelfHitScore( name      = single_hit.query_name,
                              hit_score = single_hit.hsp_score,
                              perc_query_in_align   = align_len / single_hit.query_length,
                              perc_subject_in_align = align_len / single_hit.subject_length,
                              percent_match         = single_hit.percent_match
                            )
    hit_set.append( new_score )


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
                  perc_query_in_align   = 0,
                  perc_subject_in_align = 0,
                  percent_match         = 0,
                  hit_score = 0 
                  ):
        self._name       = name
        self._other_name = other_name
        self._hit_score  = hit_score
        self._perc_q_in_align = perc_query_in_align
        self._perc_s_in_align = perc_subject_in_align
        self._perc_match      = percent_match

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
                  perc_query_in_align   = 0,
                  perc_subject_in_align = 0,
                  percent_match         = 0,
                  hit_score = 0
                ):
        super().__init__( name = name, other_name = name,
                          hit_score = hit_score,
                          perc_query_in_align = perc_query_in_align,
                          perc_subject_in_align = perc_subject_in_align,
                          percent_match = percent_match
                        )

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
                     hsp_score        = first_hsp.bits,
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

