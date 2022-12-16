#!/usr/bin/env python3

# By Jason Ladner

# Reads a set of sequences in fasta format and identifies reads that are redundant based on kmers shared with the rest of the sequences

import argparse, os, random
import numpy as np
import inout as io               #Available at https://github.com/jtladner/Modules
import fastatools as ft          #Available at https://github.com/jtladner/Modules
import kmertools as kt          #Available at https://github.com/jtladner/Modules
from collections import defaultdict


def main():

    #To parse command line
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument('-k', '--kmers', default="9,10,11,12", help='Comma-separated list of kmer sizes to consider')
    p.add_argument('-r', '--redundThresh', type=float, default=1, help='Proportion of redundant kmers to be flagged.')
    p.add_argument('-m', '--metadata', help='Tab-delimited file that can be used to link sequence names to group.')
    p.add_argument('--nameHeader', default="PeptideName", help='Header for column of metadata file containing sequence name.')
    p.add_argument('--groupHeader', default="Group", help='Header for column of metadata file containing Group info.')
    p.add_argument('--groupInfoOut', help='Optional, tsv file with info about group redundancy.')
    p.add_argument('--redRepOut', help='Optional, reduced representation output fasta.')
    p.add_argument('--rro_kRange', default="9,30", help='Range of kmer sizes to consider for reduced representation output.')
    p.add_argument('--rro_prop', default=1, type=int, help='Proportion of redundant kmers needed to be eliminated in reduced representation output.')

    reqArgs = p.add_argument_group('required arguments')
    reqArgs.add_argument('-i', '--inp', help='Either a path to directory containing designs from SW_SC.py OR a path to a file containing multiple directory paths with designs.', required=True)
#    reqArgs.add_argument('-o', '--out', help='Name for output file that will contain selected peptides from each cluster.', required=True)

    args = p.parse_args()
    
    # Parse kmer sizes
    kSizes = [int(x) for x in args.kmers.split(",")]
    
    # Read in metadata, if provided
    if args.metadata:
        metaD = io.fileDictHeader(args.metadata, args.nameHeader, args.groupHeader)
        groupD = defaultdict(list)
        totalD = defaultdict(int)
    
    # Read in sequence fasta
    fD = ft.read_fasta_dict_upper(args.inp)
    
    #Dict to hold redundant peptides
#     redunD = defaultdict(list)
#     
#     # Step through each kmer size
#     for k in kSizes:
#         
#         #Generate a dictionary of all unique kmers and their counts
#         kD = kt.kmerDictCountFasta(args.inp, k)
#         
#         #Step through each sequence individually
#         for n,s in fD.items():
#             #Get all unique kmers for this specific sequence
#             kSet = kt.kmerSet(s, k)
#             #Calculate number of redundant kmers
#             numRedund = sum([1 for kmer in kSet if kD[kmer]>1])
#             propRedund = numRedund/len(kSet)
#             # Check to see if this sequence should be flagged
#             if propRedund >= args.redundThresh:
#                 redunD[k].append(n)
#                 if args.metadata:
#                     groupD[metaD[n]].append(n)
#             totalD[metaD[n]]+=1
#     
#     # Summarize results
#     for k, v in redunD.items():
#         print(f"{k}\t{len(v)}")
#         if args.groupInfoOut and args.metadata:
#             with open(f"{k}_{args.groupInfoOut}", "w") as fout:
#                 fout.write("Group\tPropRedund\tRedundSeqs\tTotalSeqs\n")
#                 for g, l in groupD.items():
#                     fout.write(f"{g}\t{len(l)/totalD[g]}\t{len(l)}\t{totalD[g]}")
    
    # Generate set of sequences that minimize redundancy
    if args.redRepOut:
        minK, maxK = [int(x) for x in args.rro_kRange.split(",")]
        kRange = list(range(minK, maxK+1, 1))[::-1]
        rmvCount=0
        for k in kRange:
            print(f"k={k}")

            #Dictionary to hold scores
            scoreD = {}

            #Pull out counts for all kmers
            kD = kt.kmerDictCount(list(fD.values()), k)
                
            #Dictionary to hold individual kmer sets
            indivKmerSets = {}
            
            #Step through each sequence individually
            for n,s in fD.items():
                #Get all unique kmers for this specific sequence
                kSet = kt.kmerSet(s, k)
                indivKmerSets[n] = kSet
                
                #Calculate number of redundant kmers
                numRedund = sum([1 for kmer in kSet if kD[kmer]>1])
                propRedund = numRedund/len(kSet)
                if propRedund >= args.rro_prop:
                    scoreD[n] = sum([kD[each] for each in kSet])

            while len(scoreD) > 0:
                #Find the max score
                maxScore = max(scoreD.values())
                maxSeqs = [k for k,v in scoreD.items() if v==maxScore]
                toRmv = random.choice(maxSeqs)

                #Remove selected sequence from fasta dict and score dict
                del(fD[toRmv])
                del(scoreD[toRmv])

                #Adjust kD to remove counts for kmers in removed sequence
                for rK in indivKmerSets[toRmv]:
                    kD[rK]-=1
                
                #Adjust scoreD to adjust scores and remove sequences that are no longer fully redundent
                for n in list(scoreD.keys()):
                    #Calculate number of redundant kmers
                    numRedund = sum([1 for kmer in indivKmerSets[n] if kD[kmer]>1])
                    propRedund = numRedund/len(indivKmerSets[n])
                    if propRedund < args.rro_prop:
                        del(scoreD[n])
                    else:
                        scoreD[n] = sum([kD[each] for each in indivKmerSets[n]])

                rmvCount+=1
                if rmvCount%25==0:
                    print(rmvCount, len(fD))
        
        ft.write_fasta_dict(fD, args.redRepOut)
            
    # **** Write set cover style algorithm to remove sequences one by one in a way that attempts to remove those that are most redundant first
    # Start with very large kmer sizes, work down to smaller kmer sizes
    # Calculate a score for each sequence that will be the sum of the counts associated with their kmers
    # If a certain proportion of the kmers in a seq are redundant and it has the highest score among redundat seqs, it goes and scores are recalculated
    # Keep going until there aren't any full redundent peptide, switch to next k
    
###-----------------End of main()--------------------------->>>

###------------->>>

if __name__ == "__main__":
    main()
