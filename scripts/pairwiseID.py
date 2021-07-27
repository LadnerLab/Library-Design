#!/usr/bin/env python

import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules
import itertools as it
import numpy as np

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files (aligned).", nargs="+")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-o", "--out", help="Output file name. ", required=True )

    args = parser.parse_args()
    
    #Create set of characters to exclude
    exSet = set(args.exclude)

    with open(args.out, "w") as fout:
        fout.write("File\tAvgIdentity\tMinIdentity\tMaxIdentity\tAvgDivergence\tMinDivergence\tMaxDivergence\n")
        
        #Step through input files
        for eachF in args.inputs:
            
            fNames, fSeqs = ft.read_fasta_lists(eachF)
            
            ids = []
            
            #Step through each sequence pair
            for s1,s2 in it.combinations(fSeqs, 2):
                ids.append(compSeqs(s1, s2, exSet))
            
            avg = np.mean(ids)
            mn = min(ids)
            mx = max(ids)
            
            fout.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (eachF, avg, mn, mx, 1-avg, 1-mn, 1-mx))

#----------------------End of main()

def compSeqs(s1, s2, exSet):
    total=0
    same=0
    for a,b in zip(s1,s2):
        if a not in exSet and b not in exSet:
            total+=1
            if a.upper() == b.upper():
                same+=1
    return same/total

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

