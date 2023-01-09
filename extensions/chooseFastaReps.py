#!/usr/bin/env python

import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files. Facilitates batch processing. Output names will be generated using each input file name.", nargs="*")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded. By default this will be done for both the SW and SC portions of the design. However, the behavior for C residues will be different in the SW portion, when used in combination with '--swCtoS'.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-k", "--kMerSize", type=int, help="Size of kmersto use for choosing a representative.", required=True)
    reqArgs.add_argument("-o", "--out", help="Name for output fasta file.", required=True)

    args = parser.parse_args()
            
    #Create set of characters to exclude
    args.exSet = set(args.exclude)
    
    #Dictionary to hold representative sequences
    outFD = {}
    
    #Choose representatives for each input fasta
    for each in args.inputs:   #Step through each input file
        #Run the design
        n,s = chooseRep(each, args)
        outFD[n] = s
                
    # Generate concatenated output files
    ft.write_fasta_dict(outFD, args.out)
        


#----------------------End of main()

def chooseRep(inp, args):

    # Generate dict with xmer counts
    xcD = {}
    
    # Read in target sequences
    tN, tS = ft.read_fasta_lists(inp)
    
    # Read in all target Xmers
    for s in tS:
        xL = kt.kmerList(s, args.kMerSize)
        for x in xL:
            if len(set(x).intersection(args.exSet)) == 0:
                xcD[x] = xcD.get(x, 0) + 1
    
    # Score each target sequence by summing contained xmer scores. This is to choose the representative for the sliding window portion of the design
    maxScore = -1
    repS = ""
    repN = ""
    for i,s in enumerate(tS):   # Stepping through each target sequence
        theseXs = kt.kmerList(s, args.kMerSize)
        thisScore = sum([xcD[x] for x in theseXs if x in xcD])
        if thisScore > maxScore:
            maxScore = thisScore
            repS = s
            repN = tN[i]

    return repN, repS



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

