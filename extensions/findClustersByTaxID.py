#!/usr/bin/env python

import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules
import glob,os, shutil

from collections import defaultdict

def map_gbacc2tid(gbD, a2tidFL):
    gbtidD = {}
    percTarget = 10

    for file in a2tidFL:
        with open(file, "r") as fin:
            next(fin)
            for line in fin:
                cols = line.rstrip("\n").split("\t")
                if cols[0] in gbD:
                    gbtidD[cols[0]] = cols[1]
                    if len(gbtidD) == len(gbD):
                        return gbtidD
                    elif len(gbtidD)/len(gbD)*100 >= percTarget:
                        print(f"{len(gbtidD)/len(gbD)*100:.2f}% complete (percTarget={percTarget})")
                        percTarget += 10

    print(f"Did not find {len(gbD)-len(gbtidD)} accessions.")
    #Fill in missing accessions with blank values
    for k in gbD:
        if k not in gbtidD:
            gbtidD[k] = ""
            
    return gbtidD

def make_new_dir(dirpath):
    if os.path.exists(dirpath):
        if not os.path.isdir(dirpath):
            print(f"{dirpath} exists, but it is NOT a directory. Execution will be terminated.")
            return False
    else:
        os.mkdir(dirpath)
    
    return True


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--propThresh", default=0.5, type=float, help="Minimum proportion of seqs in cluster with target TaxID for cluster to be selected.")
    parser.add_argument("--copyClusters", default=False, action="store_true", help="Use this option if you want the identified clusters to be moved to a different location.")
    parser.add_argument("-w", "--warnings", default="warnings.txt", help="File to which warnings will be written.")

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument("-i", "--input", help="Tax ID(s) of interest. Tab-delimited file with one Tax ID per line. Each line should have two columns, the first is the Tax ID, the second is the base directory for the clusters for the TaxID.", required=True)
    reqArgs.add_argument("-a", "--acc2tid", help="Path to prot.accession2taxid file(s) downloaded from GenBank. If there are multiple files, just include portion of file name that is common to all", required=True)
    reqArgs.add_argument("-t", "--tidLineage", help="Path to taxidlineage.dmp file downloaded from GenBank.", required=True)
    reqArgs.add_argument("-o", "--out", help="Name of output directory, which will be created, if it doesn't already exist. If multiple Tax IDs are provided, subdirectories will be generated within this directory.", required=True)

    args = parser.parse_args()
    
    # Open output file for wanrings
    warnOut = open(args.warnings, "w")
    
    # Generate list with accession2taxid file paths
    a2tidFL = glob.glob(f"{args.acc2tid}*")
    if len(a2tidFL) == 0:
        print("WARNING: No accession2taxid files were found in the location specified! Execution will be terminated.")
        return
    else:
        print(f"Found {len(a2tidFL)} accession2taxid files")
        
    # Check to see if output directory exists. If it doesn't create it
    if not make_new_dir(args.out):
        return None
    
    # Initiate dictionaries to hold target IDs and GenBank Accessions
    targetIDs = {}
    accD = {}
    tid2common = {}
    
    # Step through each target ID
    with open(args.input, "r") as fin:
        for line in fin:
            inpCols = line.rstrip("\n").split("\t")
            tid = inpCols[0]
            tid2common[tid] = inpCols[1]
            clustDirs = inpCols[2:]
            
            # Initiate sub-dictionary to hold target IDs for this group
            targetIDs[tid] = {tid: ""}
            
            # Find all downstream IDs for the ID of interest
            # Read in TaxID lineages
            with open(args.tidLineage) as fLin:
                for line in fLin:
                    cols = line.rstrip("\n").split("\t|")
                    upstreamIDs = cols[1].split()
                    if tid in upstreamIDs:
                        targetIDs[tid][cols[0]] = ""

            print(f"{tid} ({inpCols[1]}): {len(clustDirs)} directories, {len(targetIDs[tid])} possible IDs")
    
            # Get GenBank accessions for each cluster
            accD[tid] = defaultdict(dict)
            
            fpL = []
            for clustDir in clustDirs:
                fpL += glob.glob(f"{clustDir}/*/*_id_70_*")
                
            # Step through each cluster file
            for fp in fpL:
                protName = fp.split("/")[-2]
                
                #Read cluster fasta file
                fD = ft.read_fasta_dict_upper(fp)

                # Step through each sequence name
                accL = []
                for sn in fD:
                    if "|" not in sn:
                        acc = sn.split()[0]
                        accL.append(acc)

                #Add accession list to dictionary, if any GenBank Accessions were found
                if accL:
                    accD[tid][protName][fp] = accL
                else:
                    warnOut.write(f"No GenBank accessions found for {fp}\n")
                    warnOut.write(f"{fD.keys()}\n\n")

    # Now that all accessions and IDs of interest have been identified, map accessions to TaxIDs and find clusters of interest

    # Generate dictionary with a key for every accession of interest
    allAccD = {}
    for tid, subD in accD.items():
        for pn, cf in subD.items():
            for fp, accL in cf.items():
                for ac in accL:
                    allAccD[ac] = ""
    print(f"Starting to look for {len(allAccD)} accessions")
    # Extract TaxIDs for each accession
    allAcc2tidD = map_gbacc2tid(allAccD, a2tidFL)
    print(f"Finished looking for accessions.")

    for tid, subD in accD.items():

        # Identify clusters of interest
        clustersByProt = {}

        for pn, cf in subD.items():
            clustersByProt[pn] = {}
            for fp, accL in cf.items():
                thisProp = sum([1 for x in accL if allAcc2tidD[x] in targetIDs[tid]])/len(accL)
                if thisProp >= args.propThresh:
                    clustersByProt[pn][fp] = (len(accL), thisProp)

        # Make directory for this TaxID, if it doesn't exist
        if not make_new_dir(f"{args.out}/{tid}_{tid2common[tid]}"):
            return None
    
        # Open output file for writing summary stats
        with open(f"{args.out}/{tid}_{tid2common[tid]}/{tid}_clusterSummary.tsv", "w") as fout:
            fout.write("Protein\tCluster\tNumSeqs\tNumAccessions\tPropMatched\n")
        
            # Step through each protein
            for k,v in clustersByProt.items():
                if len(v) == 0:
                    print(f"No clusters were found for {tid}:{k}")
                else:
                    # Make directory for this protein, if it doesn't exist
                    if not make_new_dir(f"{args.out}/{tid}_{tid2common[tid]}/{k}"):
                        return None

                    fD_out = {}
                    for each, info in v.items():
                        fD = ft.read_fasta_dict_upper(each)
                        fout.write(f"{k}\t{os.path.basename(each)}\t{len(fD)}\t{info[0]}\t{info[1]}\n")
                        for n, s in fD.items():
                            fD_out[n] = s

                        # Copy clusters, if requested
                        if args.copyClusters:
                            shutil.copy(each, f"{args.out}/{tid}_{tid2common[tid]}/{k}/{os.path.basename(each)}")
    
    #Close warnings file
    warnOut.close()

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

