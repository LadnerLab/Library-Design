#!/usr/bin/env python

import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules
import numpy as np
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("inputs", help="One or more input target fasta files (unaligned).", nargs="+")
    parser.add_argument("-e", "--exclude", default="X-", help="Any Xmers or yMers containing these chaarcters will be excluded.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument( '-k', '--kmer_size', help = "Comma-delimited list of kmer sizes to use for comparing sequences.", required=True )
    reqArgs.add_argument("-o", "--out", help="Output file name. ", required=True )

    args = parser.parse_args()
    
    #Create set of characters to exclude
    exSet = set(args.exclude)

    #Parse kmers
    kmers = [int(k) for k in args.kmer_size.split(",")]
    
    with open(args.out, "w") as fout:
        fout.write("File\t%s\t%s\n" % ("\t".join(["Avg%dmerProp" % k for k in kmers]), "\t".join(["Avg%dmers" % k for k in kmers])))
        
        #Step through input files
        for eachF in args.inputs:
            
            fNames, fSeqs = ft.read_fasta_lists(eachF)
            
            avgProps = []
            
            #Step through each kmer size
            for k in kmers:
                cD = kt.kmerDictCountFasta(eachF,k,filter=exSet)
                avgProps.append(np.mean(list(cD.values())))
        
    
            fout.write("%s\t%s\t%s\n" % (eachF, "\t".join(["%.3f" % (ap/len(fNames)) for ap in avgProps]), "\t".join(["%.3f" % (ap) for ap in avgProps])))
    

#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

