#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import fastatools as ft
import inout as io
import argparse, os
from collections import defaultdict

def checkIfGood(row, args):
    #Calculate and evaluate "subject ratio"
    sr = (row["Subject End"]-row["Subject Start"]+1)/row["Subject Length"]
    if sr >= args.subject_ratio:
        return True
    
    #Calculate and evaluate "query ratio"
    qr = (row["Query End"]-row["Query Start"]+1)/row["Query Length"]
    if qr >= args.query_ratio:
        return True
    
    #Determine whether the hit starts very close to the beginning of the query
    if row["Query Start"] <= args.dist_to_end and row["Alignment Length"] >= args.min_len:
        return True

    #Determine whether the hit ends very close to the end of the query
    if (row["Query Length"] - row["Query End"] + 1) <= args.dist_to_end and row["Alignment Length"] >= args.min_len:
        return True


parser = argparse.ArgumentParser()
#parser.add_argument('-u','--unassigned_filename',required=True, help="FASTA file of sequences with unassigned regions.")
parser.add_argument('-m','--metadata',required=True, help="txt file of metadata, see documentation.")
parser.add_argument('-r', '--subject_ratio', type=float, default=0.9, help='threshold for subject ratio [0.9]')
parser.add_argument('-q', '--query_ratio', type=float, default=0.9, help='threshold for query ratio [0.9]')
parser.add_argument('-d', '--dist_to_end', type=int, default=7, help='threshold for subject ratio [7]')
parser.add_argument('-l', '--min_len', type=int, default=50, help='Minimum length of hit, used in combination with --dist_to_end [50]')
parser.add_argument('--outputSubDirs', default=False, action="store_true", help='Use this flag if you want output files placed into subdirectories based on protein names.')
args=parser.parse_args()


# Read in metadata that links a blast result to a query file and a protein name
blasts=io.fileDictHeader(args.metadata,"Blast Result", "Query Fasta",)    
proteins=io.fileDictHeader(args.metadata,"Blast Result", "Protein")

# Create subdirectories, if needed
if args.outputSubDirs:
    for p in set(proteins.values()):
        if not os.path.isdir(p):
            os.mkdir(p)

#Create dictionary to keep track of clusters with no matches for certain proteins
noMatches = defaultdict(list)


# Step through each blast result file
for blast, query in blasts.items():
    
    #Read in query fasta file
    fD=ft.read_fasta_dict_upper(query)
    
    #Read blast results into a pandas dataframe
    df=pd.read_csv(blast, sep='\t',header=0)
    

    #Dictionary to hold output substrings
    outD = {}
    
    #Step through each hit and check to see if it is "good"
    for i,row in df.iterrows():
        if checkIfGood(row, args):
            #Generate name
            n = ("%s %d-%d" % (row["Query Name"], row["Query Start"], row["Query End"]))
            #Generate subsequence
            s = fD[row["Query Name"]][row["Query Start"]-1:row["Query End"]]
            #Add to output dictionary
            outD[n] = s
    
    
    #Write out fasta file with subset sequences, but only if some good hits were found
    if len(outD) > 0:
        if args.outputSubDirs:
            ft.write_fasta_dict(outD,'%s/%s_%s'%(proteins[blast], proteins[blast], os.path.basename(query)))
        else:
            ft.write_fasta_dict(outD,'%s_%s'%(proteins[blast], os.path.basename(query)))
    else:
        noMatches[os.path.basename(query)].append(proteins[blast])

if len(noMatches) >0:
    print("No matches found for the following:")
    for k,v in noMatches.items():
        print("%s:\t%s" % (k, "\t".join(v)))
