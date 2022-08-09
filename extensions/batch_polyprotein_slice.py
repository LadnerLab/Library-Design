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

def checkOverlap(hits):
    top = hits[0]
    coveredSubPos = set(range(top["Subject Start"], top["Subject End"]+1))
    other = hits[1:]
    
    toKeep = [top]
    for each in other:
        theseSubPos = set(range(each["Subject Start"], each["Subject End"]+1))
        ovlp = coveredSubPos.intersection(theseSubPos)
        if len(ovlp)/len(coveredSubPos) < 0.05:
            toKeep.append(each)
            coveredSubPos = coveredSubPos.union(theseSubPos)
    
    return(toKeep)

# Checks to make sure the relative orders of hits in the query and subject are the same
# If they aren't, it will remove hits from the hit (lower scoring) until condordance is reached
def checkConcordance(hits):
    qStarts = [(h["Query Start"], i) for i, h in enumerate(hits)]
    sStarts = [(h["Subject Start"], i) for i, h in enumerate(hits)]

    while [x[1] for x in sorted(qStarts)] != [x[1] for x in sorted(sStarts)]:
        hits = hits[:-1]
        del(qStarts[-1])
        del(sStarts[-1])

    return hits

parser = argparse.ArgumentParser()
#parser.add_argument('-u','--unassigned_filename',required=True, help="FASTA file of sequences with unassigned regions.")
parser.add_argument('-m','--metadata',required=True, help="txt file of metadata, see documentation.")
parser.add_argument('-e', '--expect', type=float, default=0.000001, help='threshold for Hsp Expect [0.000001]')
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
    
    #Remove any hits that don't meet the "Hsp Expect" threshold
    df = df[df["Hsp Expect"] <= args.expect]
    
    ###Combine Hsps for the same query/subject pair
    
    #Count the number of times each query subject pair appears in the output
    qsPairs = defaultdict(int)
    for i,row in df.iterrows():
        qsPairs[(row["Query Name"], row["Subject Name"])] += 1
#         if qsPairs[(row["Query Name"], row["Subject Name"])] != 1:
#             print ((row["Query Name"], row["Subject Name"]), qsPairs[(row["Query Name"], row["Subject Name"])])
    
    #Generate new combo dataframe, with a subset of the columns
    #dfComboD = {k:[] for k in ["Query Name", "Query Start", "Query End", "Query Length", "Subject Name", "Subject Start", "Subject End", "Subject Length", "Alignment Length"]}
    toCombine = defaultdict(list)
    dfCombo = pd.DataFrame()
    
    # Step through each hit and either write to dfCombo or save info in toCombine for potential combination
    for i,row in df.iterrows():
        if qsPairs[(row["Query Name"], row["Subject Name"])] == 1:
            dfCombo = dfCombo.append(row)
        else:
            toCombine[(row["Query Name"], row["Subject Name"])].append(row)
            
    
    
    for k,v in toCombine.items():
        tempD = {}

        #Check to make sure hits to combine don't have larger overlaps of the same portion of the reference
        v = checkOverlap(v)
        
        #Check to make sure these is concordance among the good hits (i.e., relative order in subject and query are the same)
        v = checkConcordance(v)
        
        #If additional hits all overlap substantially with the top hit
        if len(v) == 1:
            dfCombo = dfCombo.append(v[0])

        #If there are multiple hits that are non-overlapping
        else:
        
            #Categories that can just be take from the first Hsp
            for each in ["Query Name", "Subject Name", "Query Length", "Subject Length"]:
                tempD[each] = v[0][each]

            #Categories for which we'll take the minimum
            for each in ["Query Start", "Subject Start"]:
                tempD[each] = min([r[each] for r in v])

            #Categories for which we'll take the maximum
            for each in ["Query End", "Subject End"]:
                tempD[each] = max([r[each] for r in v])
        
            #Categories that need to be calculated
            tempD["Alignment Length"] = tempD["Query End"] - tempD["Query Start"] + 1
        
            #Add new combo entry to dataframe
            dfCombo = dfCombo.append(tempD, ignore_index=True)
        
            #This part is temporary to check that things are working
#             for each in v:
#                 print(each)
#             print("")
#             print(tempD)
#             print("\n")

#     print(dfCombo.shape)
    
    #Step through each hit and check to see if it is "good"
    for i,row in dfCombo.iterrows():
        if checkIfGood(row, args):
            #Generate name
            n = ("%s %d-%d" % (row["Query Name"], int(row["Query Start"]), int(row["Query End"])))
            #Generate subsequence
            s = fD[row["Query Name"]][int(row["Query Start"])-1:int(row["Query End"])]
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
