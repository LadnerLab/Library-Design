#!/usr/bin/env python3

import fastatools as ft  # Available here: https://github.com/jtladner/Modules
import inout as io  # Available here: https://github.com/jtladner/Modules
import argparse
from collections import defaultdict


def main():
    arg_parser = argparse.ArgumentParser( description = "Extract sequences from a GenBank fasta file based on user provided keywords.", formatter_class=argparse.ArgumentDefaultsHelpFormatter )

    arg_parser.add_argument('--accHeader', default="Accession", help = "Header name for accession number in 'info' file.")
    arg_parser.add_argument('--typeHeader', default="SubType", help = "Header name for subtype in 'info' file.")
    arg_parser.add_argument('--nameHeader', default="SubName", help = "Header name for subname in 'info' file.")
    arg_parser.add_argument('--orgHeader', default="Organism", help = "Header name for subname in 'info' file.")
    arg_parser.add_argument('--cats', default="host", help = "Categories to parse for keywords. Comma-separated")
    arg_parser.add_argument('--delim', default="|", help = "Delimiter for parsing subtype and subname.")
    arg_parser.add_argument('--delimList', default=[" ", "_", "(", ")", "|", ","], help = "Delimiters for parsing host and organism strings.")

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-f', '--fasta', help = "Fasta file contianing the protein sequence(s) downloaded from GenBank. First piece of info in the name should be the accession", required=True )
    reqArgs.add_argument( '-i', '--info', help = "Tab-delimited file containing a column with accession and two additional columns: SubType and SubName.", required=True )
    reqArgs.add_argument( '-k', '--keys', help = "File containing list of keywords of interest.", required=True )
    reqArgs.add_argument( '-e', '--excl', help = "File containing list of keywords that should be excluded.", required=True )
    reqArgs.add_argument( '-o', '--output', help = "Name for output fasta file containing fastas of interest.", required=True )

    args = arg_parser.parse_args()
    
    # Prep categories
    catD = {k:"" for k in args.cats.split(",")}
    
    # Read in data from info file
    subTypeD = io.fileDictHeader(args.info, args.accHeader, args.typeHeader)
    subNameD = io.fileDictHeader(args.info, args.accHeader, args.nameHeader)
    orgD = io.fileDictHeader(args.info, args.accHeader, args.orgHeader)
    
    # Read in seqs from fasta
    fD = ft.read_fasta_dict_upper(args.fasta)
    
    # Make dict linking accession to fasta seq name
    acc2fasD = {n.split()[0]:n for n in fD}
    
    # Read in keywords
    kwD = io.fileEmptyDict(args.keys, header=False)
    kwD = {k.upper():v for k,v in kwD.items()}

    exD = io.fileEmptyDict(args.excl, header=False)
    exD = {k.upper():v for k,v in exD.items()}
    
    # Set up containers for outputs
    outD = {}
    excludeD = defaultdict(int)
    
    # Step through each seq and look for keywords
    for acc, sn in acc2fasD.items():
        match = 0
        exclMatch = 0
        
        if acc in orgD:
            org = orgD[acc]
            infoTypes = {k:i for i,k in enumerate(subTypeD[acc].split(args.delim))}
            infoNames = subNameD[acc].split(args.delim)
        
            #Parse host info, etc. 
            for k,i in infoTypes.items():
                if k in catD:
                    thisInfo = [x.upper() for x in splitMultiDelim(infoNames[i], args.delimList)]
                    for each in thisInfo:
                        if each in kwD:
                            match += 1
                        if each in exD:
                            exclMatch += 1

                    if exclMatch == 0 and match ==0:
                        excludeD[infoNames[i]] += 1

            # Parse organism info
            thisInfo = [x.upper() for x in splitMultiDelim(org, args.delimList)]
            for each in thisInfo:
                if each in kwD:
                    match += 1
                if each in exD:
                    exclMatch += 1

        else:
            thisInfo = [x.upper() for x in splitMultiDelim(sn, args.delimList)]
            for each in thisInfo:
                if each in kwD:
                    match += 1
                if each in exD:
                    exclMatch += 1

        if match > exclMatch:
            seqName = acc2fasD[acc]
            outD[seqName] = fD[seqName]
                    
                    
    # Write out fasta file with seqs of interest
    ft.write_fasta_dict(outD, args.output)
    
    # Print to the screen the info for seqs that were excluded
    for k,v in excludeD.items():
        print(k,v)

def splitMultiDelim(string, delims):

    if len(delims) > 1:
        for each in delims[1:]:
            string = string.replace(each, delims[0])

    spl = string.split(delims[0])
    return [x for x in spl if x]

if __name__ == '__main__':
    main()
