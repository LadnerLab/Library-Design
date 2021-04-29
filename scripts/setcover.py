#!/usr/bin/env python

import argparse, random
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script uses a set cover algorithm to design peptides that optimize linear epitope coverage
# This is analogous to the design process implemented in C++ here: https://github.com/LadnerLab/C-KmerOligo

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--pre", help="Comma-sep list of fasta files containing previously designed peptides. Xmers contained in these sequecnes will not contribute to Ymer scoring in design.")
    parser.add_argument("-t", "--target", default=1, type=float, help="Target ymer coverage. Algorithm will continue until at least this proportion of total Ymers are in the design. If '--pre' option is used, Ymers in these predesigned peptides will also be considered in this threshold.")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")
#    parser.add_argument("--includeTerminalDashes", default=True, action="store_false", help="By default, terminal '-' characters will not be considered in consensus generation.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-i", "--inp", help="Input file name. Should contain target protein sequences from which to design peptides.", required=True)
    reqArgs.add_argument("-o", "--out", help="Output file name. Will be a list of peptides, 1 per line.", required=True)
    reqArgs.add_argument("-x", "--xMerSize", type=int, help="Size of Xmers, which represent potential linear epitopes contained within peptides/Yemrs.", required=True)
    reqArgs.add_argument("-y", "--yMerSize", type=int, help="Size of Ymers, which represent potential linear epitopes contained within peptides/Yemrs.", required=True)

    args = parser.parse_args()
    
    #Create set of characters to exclude
    exSet = set(args.exclude)
    
    # Generate dict with kmer counts
    ycD = defaultdict(int)
    
    tN, tS = ft.read_fasta_lists(args.inp)
    for s in tS:
        yL = kt.kmerList(s, args.yMerSize)
        for y in yL:
            if len(set(y).intersection(exSet)) == 0:
                ycD[y]+=1
    
    #Save count of total ymers in targets
    totalY = len(ycD)
    
    # If pre-designed peptides are provided, remove any contained ymers from the ycD
    if args.pre:
        for each in args.pre.split(","):
            pN, pS = ft.read_fasta_lists(each)
            for s in pS:
                yL = kt.kmerList(s, args.yMerSize)
                for y in yL:
                    if y in ycD:
                        del(ycD[y])
            
    # Read in all xMers in targets
    xsD = {}
    for s in tS:
        xL = kt.kmerList(s, args.xMerSize)
        for x in xL:
            if len(set(x).intersection(exSet)) == 0:
                xsD[x] = 0
    
    # Design peptides
    newPeps = []
    
    while (1-(len(ycD)/totalY)) < args.target:
        
        thisPep = choosePep(xsD, ycD, args)
        newPeps.append(thisPep)
        
        #Remove selected peptide from xsD
        del(xsD[thisPep])
        
        #Remove covered yMers from ycD
        for eachY in kt.kmerList(thisPep, args.yMerSize):
            if eachY in ycD:
                del(ycD[eachY])
        
    with open(args.out, "w") as fout:
        fout.write("\n".join(newPeps))

#----------------------End of main()

def choosePep(xsD, ycD, args):
    #Calculate scores for xMers
    for x in xsD:
        theseYs = kt.kmerList(x, args.yMerSize)
        xsD[x] = sum([ycD[y] for y in theseYs if y in ycD])

    #Dict by score
    scoreD = defaultdict(list)
    for k,v in xsD.items():
        scoreD[v].append(k)
    
    #Choose peptide
    thisMax = max(scoreD.keys())
    
    thisChoice = random.choice(scoreD[thisMax])
    
    return thisChoice


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

