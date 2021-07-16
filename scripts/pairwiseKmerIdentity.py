#!/usr/bin/env python

import argparse
import itertools as it
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files (unaligned).", nargs="+")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument( '-k', '--kmer_size', help = "kmer size to use for comparing sequences.", default = 9, type = int, required=True )
    reqArgs.add_argument("-o", "--out", help="Output file name. ", required=True )

    args = parser.parse_args()
    
    #Create set of characters to exclude
    exSet = set(args.exclude)
    
    with open(args.out, "w") as fout:
        fout.write("File\tAvgPropShared\tMinPropShared\tMaxPropShared\n")
        
        #Step through input files
        for eachF in args.inputs:
            
            #Read in seqs in file
            names, seqs = ft.read_fasta_lists(eachF)
            
            propIDs = []
            
            for s1, s2 in it.combinations(seqs, 2):
                propIDs.append(kt.compSeqs(s1, s2, args.kmer_size, filter=exSet))
        
            fout.write("%s\t%.3f\t%.3f\t%.3f\n" % (eachF, sum(propIDs)/len(propIDs), min(propIDs), max(propIDs)))
    

#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

