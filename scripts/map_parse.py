#!/usr/bin/env python3
import argparse
import sys, re

import protein_oligo_library as oligo

def main():
    arg_parser = argparse.ArgumentParser( description = "Script to parse information given by validate_design's epitope map output" )

    arg_parser.add_argument( '-m', '--map', help = "Name of file containing epitope map" )
    arg_parser.add_argument( '-o', '--output', default = "parsed_map",
                             help = (
                                      "Identifier to preface output files."
                                      " Output will be written in a tab-delimited format"
                                    )
                           )
    arg_parser.add_argument( '-t', '--tax_db',
                             help = "Name of file containing mappings of taxonomic id -> rank data"
                           )
    arg_parser.add_argument( '-v', '--verbose',
                             help = "Flag to add if output should be written to STDOUT"
                           )
    arg_parser.add_argument( '-g', '--gap_file',
                             help = ( "File containing mappings of taxid->rank for "
                                      "use in filling gaps."
                                    )
                           )
    arg_parser.add_argument( '-r', '--reference',
                             help = "Fasta reference dataset used to create library"
                           )

    arg_parser.add_argument( '-c', '--code',
                             help = "Tab-delimited file needed if the map file contains coded names"
                           )


    args = arg_parser.parse_args()

    oligo_centric_table    = None
    sequence_centric_table = None
    species_centric_table  = None

    map_dict   = {}
    taxid_dict = {}
    gap_dict   = {}

    missing_id_key = {
                       10969	: 444185,
                       11619	: 2169991,
                       11630	: 2169993,
                       11806	: 353765,
                       45218	: 2169996,
                       45222	: 2169994,
                       45709	: 2169992,
                       489502	: 10407,
                       587201	: 10255,
                       587202	: 10255,
                       587203	: 10255,
                       1173522	: 11723,
                       1554474	: 1511807,
                       1554482	: 1330068,
                       1554483	: 1330491,
                       1554492	: 1330066,
                       1554494	: 1307800,
                       1554498	: 1511784,
                       1559366	: 1513237,
                       1560037	: 1131483
                     }


    #Read in code, if a code file is provided
    if args.code:
        decode={}
        with open(args.code, "r") as fin:
            for line in fin:
                cols=line.strip("\n").split("\t")
                decode[cols[1]]=cols[0]

    for line in open( args.gap_file, 'r' ):
        line = line.split( '|' )
        gap_dict[ line[ 0 ] ] = line[ 1 ].strip()
    for line in open( args.tax_db, 'r' ):
        line = line.split( '|' )
        taxid_dict[ line[ 0 ].strip() ] = [ item.strip() for item in line[ 1:: ] ]

    oligo.fill_tax_gaps( taxid_dict, gap_dict )
    
    oligo_file = open( args.output + "_oligo.tsv", 'w+' )
    oligo_file.write( 
                      "Oligo Name\tNum Sequences Share 7-mer\tNum Species Share 7-mer\t"
                      "Num Genera Share 7-mer\t"
                      "Num Families Share 7-mer\t"
                      "Species Covered by 7-mer\t"
                      "Genera Covered by 7-mer\t"
                      "Families Covered by 7-mer\n"
                    )

    missing_ids = set()
    oligo_file.close()
    for line in open( args.map, 'r' ):

        oligo_file = open( args.output + "_oligo.tsv", 'a' )
        line = line.strip("\n").split( '\t' )

        if args.code:
            taxids = set( [ get_OXXids_from_name( decode[item] ) for item in line[ 1 ].split( '~' ) ] )
        else:
            taxids = set( [ get_OXXids_from_name( item ) for item in line[ 1 ].split( '~' ) ] )
        
#        taxids = oligo.create_valid_taxids( taxids, missing_id_key )

        current_entry = line[ 0 ].strip() + '\t'

        current_species = set([x[0] for x in taxids if x[0]])
        current_genus = set([x[1] for x in taxids if x[1]])
        current_family = set([x[2] for x in taxids if x[2]])

        JOIN_CHAR = '~'

        current_entry += "%d\t%d\t%d\t%d\t" % ( len( line[ 1 ].split( '~' ) ), len( current_species ),
                                            len( current_genus ),
                                            len( current_family )
                                          )
        current_entry += "%s\t" % JOIN_CHAR.join( current_species )
        current_entry += "%s\t" % JOIN_CHAR.join( current_genus   )
        current_entry += "%s"   % JOIN_CHAR.join( current_family  )
        current_entry += '\n'

        oligo_file.write( current_entry )

        oligo_file.close()

    print( "Missing ids: %s" % ",".join( missing_ids ) )
    # Generate species- genus- and family-centric files
    
    famDict = {}
    genDict = {}
    spDict = {}
    
    famDictSpec = {}
    genDictSpec = {}
    spDictSpec = {}
    
    fin = open( args.output + "_oligo.tsv", 'r' )

    linecount=0
    for line in fin:
        linecount+=1
        if linecount==1:
            headers = line.strip("\n").split("\t")
        else:
            name, numSeqs, numSp, numGen, numFam, namesSp, namesGen, namesFam = line.strip("\n").split("\t")
            #species-centric
            for each in namesSp.strip(',').split(JOIN_CHAR):
                if len( each ) > 0:
                    spDict[each] = spDict.get(each, 0) + 1
                    if len(namesSp.strip(',').split(JOIN_CHAR))==1:
                        spDictSpec[each] = spDictSpec.get(each, 0) + 1
            #genus-centric
            for each in namesGen.strip(',').split(JOIN_CHAR):
                if len( each ) > 0:
                    genDict[each] = genDict.get(each, 0) + 1
                    if len(namesGen.strip(',').split(JOIN_CHAR))==1:
                        genDictSpec[each] = genDictSpec.get(each, 0) + 1
            #family-centric
            for each in namesFam.strip(',').split(JOIN_CHAR):
                if len( each ) > 0:
                    famDict[each] = famDict.get(each, 0) + 1
                    if len(namesFam.strip(',').split(JOIN_CHAR))==1:
                        famDictSpec[each] = famDictSpec.get(each, 0) + 1
    
    print (len(spDict), len(spDictSpec))
    print (len(genDict), len(genDictSpec))
    print (len(famDict), len(famDictSpec))
    
    #Write out species-centric
    fout = open("species-centric.txt", "w")
    fout.write("Species\tTotalOligos\tSpecificOligos\n")
    for sp,count in spDict.items():
        if sp in spDictSpec: specific = spDictSpec[sp]
        else: specific = 0
        fout.write("%s\t%d\t%d\n" % (sp, count, specific))
    fout.close()
    
    #Write out genus-centric
    fout = open("genus-centric.txt", "w")
    fout.write("Genus\tTotalOligos\tSpecificOligos\n")
    for gen,count in genDict.items():
        if gen in genDictSpec: specific = genDictSpec[gen]
        else: specific = 0
        fout.write("%s\t%d\t%d\n" % (gen, count, specific))
    fout.close()
    
    #Write out family-centric
    fout = open("family-centric.txt", "w")
    fout.write("Family\tTotalOligos\tSpecificOligos\n")
    for fam,count in famDict.items():
        if fam in famDictSpec: specific = famDictSpec[fam]
        else: specific = 0
        fout.write("%s\t%d\t%d\n" % (fam, count, specific))
    fout.close()

def get_OXXids_from_name(name):
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_ids = oxpat.search(name)
    return tax_ids.group(2),tax_ids.group(3),tax_ids.group(4)

if __name__ == '__main__':
    main()
