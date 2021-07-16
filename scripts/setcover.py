#!/usr/bin/env python

import argparse, random, os
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script uses a set cover algorithm to design peptides that optimize linear epitope coverage
# This is analogous to the design process implemented in C++ here: https://github.com/LadnerLab/C-KmerOligo

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files. Facilitates batch processing. Output names will be generated using each input file name.", nargs="*")
    parser.add_argument("-i", "--inp", help="Input file name. Should contain target protein sequences from which to design peptides. Can be used along with -o if designing for a single target set.")
    parser.add_argument("-o", "--out", help="Output file name. Will be a list of peptides, 1 per line. Can be used along with -i if designing for a single target set.")
    parser.add_argument("-p", "--pre", help="Comma-sep list of fasta files containing previously designed peptides. Xmers contained in these sequecnes will not contribute to Ymer scoring in design.")
    parser.add_argument("-t", "--target", default=1, type=float, help="Target ymer coverage. Algorithm will continue until at least this proportion of total Ymers are in the design. If '--pre' option is used, Ymers in these predesigned peptides will also be considered in this threshold.")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")
    parser.add_argument("-u", "--summary", help="Name for a tab-delimited output file summarizing the number of peptides designed for each input set of targets.")
#    parser.add_argument("--outputXmerTables", default=False, action="store_true", help="Use this flag to write out Xmer tables pre- and post- removal of Xmers from pre-selected Ymers.")
#    parser.add_argument("--includeTerminalDashes", default=True, action="store_false", help="By default, terminal '-' characters will not be considered in consensus generation.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-x", "--xMerSize", type=int, help="Size of Xmers, which represent potential linear epitopes contained within peptides/Yemrs.", required=True)
    reqArgs.add_argument("-y", "--yMerSize", type=int, help="Size of Ymers, which represent potential peptides for inclusion in the assay.", required=True)

    args = parser.parse_args()
    
    #Create set of characters to exclude
    args.exSet = set(args.exclude)
    
    # Open output summary file for writing, if requested
    if args.summary:
        fout = open(args.summary, "w")
        fout.write("File\tNumPeps\n")
    
    #Run set cover analyses
    for each in args.inputs:
        numPep = design(each, "%s_SC-x%d-y%d.fasta" % (os.path.basename(each), args.xMerSize, args.yMerSize), args)

        if args.summary:
            fout.write("%s\t%d\n" % (each, numPep))
    
    if args.inp and args.out:
        numPep = design(args.inp, args.out, args)

        if args.summary:
            fout.write("%s\t%d\n" % (args.inp, numPep))

    
#----------------------End of main()

def design(inp, out, args):

    # Generate dict with xmer counts
    xcD = defaultdict(int)
    
    tN, tS = ft.read_fasta_lists(inp)
    for s in tS:
        xL = kt.kmerList(s, args.xMerSize)
        for x in xL:
            if len(set(x).intersection(args.exSet)) == 0:
                xcD[x]+=1
    
    # Write out tsv with xmer counts, if requested
#     if args.outputXmerTables:
#         writeXmerDict(xcD, "initialXmerCounts.tsv")
    
    #Save count of total xmers in targets
    totalX = len(xcD)
    
    # If pre-designed peptides are provided, remove any contained xmers from the xcD
    if args.pre:
        for each in args.pre.split(","):
            pN, pS = ft.read_fasta_lists(each)
            for s in pS:
                xL = kt.kmerList(s, args.xMerSize)
                for x in xL:
                    if x in xcD:
                        del(xcD[x])
        
        # Write out tsv with xmer counts, if requested
#         if args.outputXmerTables:
#             writeXmerDict(xcD, "preRemovedXmerCounts.tsv")

    # Read in all yMers in targets
    ysD = {}
    yNameD = {}
    for i,s in enumerate(tS):
#    for s in tS:
        yL = kt.kmerList(s, args.yMerSize)
        for j, y in enumerate(yL):
#        for y in yL:
            if len(set(y).intersection(args.exSet)) == 0:
                ysD[y] = 0
                yNameD[y] = "%s_%04d" % (tN[i], j)
    
    # Design peptides
    newPeps = []
    newNames = []
    
    while (1-(len(xcD)/totalX)) < args.target:
        
        thisPep = choosePep(ysD, xcD, args)
        thisName = yNameD[thisPep]
        newPeps.append(thisPep)
        newNames.append(thisName)
        
        #Remove selected peptide from ysD
        del(ysD[thisPep])
        
        #Remove covered xMers from xcD
        for eachX in kt.kmerList(thisPep, args.xMerSize):
            if eachX in xcD:
                del(xcD[eachX])
        
    # Write out peptides
    ft.write_fasta(newNames, newPeps, out)

    return len(newPeps)

def choosePep(ysD, xcD, args):
    #Calculate scores for xMers
    for y in ysD:
        theseXs = kt.kmerList(y, args.xMerSize)
        ysD[y] = sum([xcD[x] for x in theseXs if x in xcD])

    #Dict by score
    scoreD = defaultdict(list)
    for k,v in ysD.items():
        scoreD[v].append(k)
    
    #Choose peptide
    thisMax = max(scoreD.keys())
    
    thisChoice = random.choice(scoreD[thisMax])
    
    return thisChoice

def writeXmerDict(xD, outname):
    with open(outname, "w") as fout:
        fout.write("Xmer\tCount\n")
        for k,v in xD.items():
            fout.write("%s\t%d\n" % (k, v))

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

