#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import fastatools as ft
import inout as io
import argparse, os

parser = argparse.ArgumentParser()
#parser.add_argument('-u','--unassigned_filename',required=True, help="FASTA file of sequences with unassigned regions.")
parser.add_argument('-m','--metadata',required=True, help="txt file of metadata, see documentation.")
parser.add_argument('-r', '--subject_ratio', type=float, required=True, help='threshold for subject ratio (ex .75)')
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

# Step through each blast result file
for blast, query in blasts.items():
    
    #Read in query fasta file
    fD=ft.read_fasta_dict_upper(query)
    
    #Read blast results into a pandas dataframe
    df=pd.read_csv(blast, sep='\t',header=0)
    
    #Add column for subject ratio
    df["Subject_Ratio"]=(df["Subject End"]-df["Subject Start"]+1)/df["Subject Length"]
    #Remove any rows that don't meet the specified subject ratio
    df=df[df["Subject_Ratio"]>=args.subject_ratio]

    #Dictionary to hold output substrings
    outD = {}
    
    #Step through each good hit and add to output dictionary
    for i,row in df.iterrows():
        #Generate name
        n = ("%s %d-%d" % (row["Query Name"], row["Query Start"], row["Query End"]))
        #Generate subsequence
        s = fD[row["Query Name"]][row["Query Start"]-1:row["Query End"]]
        #Add to output dictionary
        outD[n] = s
    
    #Write out fasta file with subset sequences, but only if some good hits were found
    if len(outD) > 1:
        if args.outputSubDirs:
            ft.write_fasta_dict(outD,'%s/%s_%s'%(proteins[blast], proteins[blast], os.path.basename(query)))
        else:
            ft.write_fasta_dict(outD,'%s_%s'%(proteins[blast], os.path.basename(query)))
    else:
        print("No matches")