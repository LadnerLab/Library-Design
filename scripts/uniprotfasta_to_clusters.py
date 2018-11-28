#!/usr/bin/env python3
import argparse                       # for parsing command-line arguments
import os                             # For creating directories
import sys                            # for exiting upon failure
import re                             # Regular expression functionality for extracting TaxIDs
import random                         # To randomly subsample protein sequences
import shutil                         # To copy files
#from enum import Enum                 # for handling errors



def main():

    parser = argparse.ArgumentParser( description = "Script to convert UniprotKB/Swiss-Prot "
                                                        "entries to a fasta file. The names for "
                                                        "each entry can be specified from the tags "
                                                        "found in each entry."
                                        )

    parser.add_argument( '-i', '--ids', default = None, type = str,
                         help = "File with target TaxIDs. "
                                "One ID per line. "
                       )
    parser.add_argument( '-s', '--sprot', default = None, type = str,
                         help = "Fasta file of swissprot entries. "
                       )
    parser.add_argument( '-t', '--trembl', default = None, type = str,
                         help = "Fasta file of trembl entries. "
                       )
    parser.add_argument( '-x', '--max', default = 5000, type = int,
                         help = "Default # of proteins to include per species. "
                       )
    parser.add_argument( '-k', '--kSize', default = 9, type = int,
                         help = "Size of kmers to calculate. "
                       )
    parser.add_argument( '-m', '--maxKmers', default = 500000, type = int,
                         help = "Max # kmers allowed in a cluster. "
                       )

    parser.add_argument( '-d', '--dir', type = str,
                             help = "Name of directory to create to hold cluster fasta files. "
                           )
                           
    parser.add_argument( '-l', '--ranked_lineage', default = None, type = str,
                         help = "Map containing taxid|tax_info pairings. "
                                "Needs to be provided if you want cluster files "
                                "to be name according to the actual taxonomic name, "
                                "instead of just the ID #."
                       )

    args = parser.parse_args()

#    args_result = validate_args( args )

#    if args_result != ArgResults.NO_ERR.value:
#        report_error( args_result )
#        sys.exit( 1 )
    
    
    #Make dict to link taxIDs to names
    #Will return None, if no file has been provided
    taxmap = get_taxmap_from_file( args.ranked_lineage )
    
    #Read in taxid names of interest
    try: 
        IDict=filedict(args.ids)
    except: 
        if args.ids:    print ("%s doesn't exist!" % (args.ids))
        else: print ("You must provide a target ids file!")
        sys.exit(1)
    
    if not os.path.exists(args.dir):
        os.mkdir(args.dir)                   #Generate main directory to hold output
        os.mkdir("%s/species" % args.dir)    #Generate subdirectory for species clusters
        os.mkdir("%s/genus" % args.dir)      #Generate subdirectory for genus clusters
        os.mkdir("%s/family" % args.dir)     #Generate subdirectory for family clusters
        os.mkdir("%s/design" % args.dir)     #Generate subdirectory for clusters to be used for design
    else: 
        print ("%s already exists!" % (args.dir))
        sys.exit(1)
    
    # parse the swissprot fasta file
    sSeqDict = read_fasta_select(args.sprot, IDict)

    # parse the trembl fasta file
    tSeqDict = read_fasta_select(args.trembl, IDict)

    #Down select
    toCluster={}
    num2fill={}
    #First choose sprot representatives
    for family,genera in sSeqDict.items():
        toCluster[family]={}
        for genus,species in genera.items():
            toCluster[family][genus]={}
            for s,sd in species.items():
                toCluster[family][genus][s]={}
                if len(sd['names'])<=args.max:
                    toCluster[family][genus][s]=sd
                    num2fill[(family,genus,s)]=args.max-len(sd['names'])
                else:
                    names,seqs = choose_random(sd,args.max)
                    toCluster[family][genus][s]={"names":names, "seqs":seqs}
                    num2fill[(family,genus,s)]=0

    del(sSeqDict)      #Cleaning up

    #Then add seqs from Trembl
    for family,genera in tSeqDict.items():
        if family not in toCluster: toCluster[family]={}
        for genus,species in genera.items():
            if genus not in toCluster[family]: toCluster[family][genus]={}
            for s,sd in species.items():
                if s not in toCluster[family][genus]: 
                    toCluster[family][genus][s]={"names":[], "seqs":[]}
                    num2fill[(family,genus,s)]=args.max
                if len(sd['names'])<=num2fill[(family,genus,s)]:
                    toCluster[family][genus][s]["names"] = toCluster[family][genus][s]["names"] + sd["names"]
                    toCluster[family][genus][s]["seqs"] = toCluster[family][genus][s]["seqs"] + sd["seqs"]
                elif num2fill[(family,genus,s)]:
                    names,seqs = choose_random(sd,num2fill[(family,genus,s)])
                    toCluster[family][genus][s]["names"] = toCluster[family][genus][s]["names"] + names
                    toCluster[family][genus][s]["seqs"] = toCluster[family][genus][s]["seqs"] + seqs

    del(tSeqDict)      #Cleaning up

    print("Species\tGenus\tFamily\tSpeciesID\tGenusID\tFamilyID\t#Seqs\tSumSeqLen")
    
    fout = open("%s/design/design_cluster_info.txt" % args.dir, "w")    #To hold info on kmer #s in clusters
    fout.write("ClusterFile\t#%dmers\n" % args.kSize)
    
    #First, write family-level clusters
    tooBigFam={}
    famClustNames={}
    for family,genera in toCluster.items():
        names=[]
        seqs=[]
        for genus,species in genera.items():
            for s,sd in species.items():
                names += sd['names']
                seqs += sd['seqs']
                print("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d" % (return_taxname(s, taxmap), return_taxname(genus, taxmap), return_taxname(family, taxmap), s, genus, family, len(sd['seqs']), sum([len(x) for x in sd['seqs']])))
        if family:
            thisout="%s/family/%s.fasta" % (args.dir,return_taxname(family, taxmap))
            write_fasta(names,seqs,thisout)
            famClustNames[family]=thisout
        else: #If there is no family assignment, treat as if family is too big, this will lead this group to be clustered by genus
            tooBigFam[family]=""

    #Check whether a family cluster is too big (in terms of # of unique kmers)
    for fam,each in famClustNames.items():
        kmers = fasta_kmers(each, args.kSize)
        if len(kmers)>args.maxKmers:
            tooBigFam[fam]=""
        else:
            shutil.copy(each,"%s/design" % args.dir)
#            print (each,len(kmers))
            fout.write("%s\t%d\n" % (each.split("/")[-1], len(kmers)))

    #Then genus-level clusters
    tooBigGen={}
    genClustNames={}
    for family,genera in toCluster.items():
        for genus,species in genera.items():
            names=[]
            seqs=[]
            for s,sd in species.items():
                names += sd['names']
                seqs += sd['seqs']
            if genus:
                thisout="%s/genus/%s.fasta" % (args.dir,return_taxname(genus, taxmap))
                write_fasta(names,seqs,thisout)
            else:
                tooBigGen[genus]=""
            if family in tooBigFam:
                genClustNames[genus]=thisout

    #Check whether a genus cluster is too big (in terms of # of unique kmers)

    for gen,each in genClustNames.items():
        kmers = fasta_kmers(each, args.kSize)
        if len(kmers)>args.maxKmers:
            tooBigGen[gen]=""
        else:
            shutil.copy(each,"%s/design" % args.dir)
#            print (each,len(kmers))
            fout.write("%s\t%d\n" % (each.split("/")[-1], len(kmers)))
            

    #Then species-level clusters
    spClustNames={}
    for family,genera in toCluster.items():
        for genus,species in genera.items():
            for s,sd in species.items():
                thisout="%s/species/%s.fasta" % (args.dir,"_".join(return_taxname(s, taxmap).split()))
                write_fasta(sd['names'],sd['seqs'],thisout)
                if genus in tooBigGen:
                    spClustNames[s]=thisout

#    sys.exit( 1 )

##########------------------------------------------------->>>>>>>>>

def fasta_kmers(fasta, k):
    seq_count=0
    names, seqs = read_fasta_lists(fasta)
    kmer_set = set()
    del(names)
    for s in seqs:
        seq_count+=1
        for i in range(0,len(s)-k+1,1):
            this_kmer = s[i:i+k].upper()
            if 'X' not in this_kmer:
                kmer_set.add(this_kmer)
    return kmer_set

def return_taxname(id, dict):
    if id in dict:
        return dict[id]
    else:
        return id

def choose_random(seqDict,num):
    rand_ind = random.sample(list(range(len(seqDict['names']))), num)
    return [seqDict['names'][i] for i in rand_ind], [seqDict['seqs'][i] for i in rand_ind]

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_select(file, IDict):
    count=0
    #oxpat = re.compile("OX=(\d+)")         #Initialize pattern for parsing OX TaxID from name
    oxxpat = re.compile("OXX=(\d+),(\d*),(\d*),(\d*)")         #Initialize pattern for parsing OXX TaxID from name

    seq_dict={}

    include=0
    with open(file, 'r') as fin:
        seq=''
        for line in fin:
            line=line.strip()
            if line and line[0] == '>':                #indicates the name of the sequence

                if count>1 and include:
                    # Only grab the first item in the line,
                    # protects against SEQ followed by anything else
                    seq_dict[fam][gen][sp]["seqs"].append(seq.strip().split()[0])

                count+=1
                n=line[1:]   #Sequence name
                oxx = oxxpat.search(n)
                if oxx:
                    id,sp,gen,fam = oxx.group(1),oxx.group(2),oxx.group(3),oxx.group(4)
                    if sum([1 for x in [id,sp,gen,fam] if x in IDict]):
                        include=1
                        if fam not in seq_dict:
                            seq_dict[fam]={}
                        if gen not in seq_dict[fam]:
                            seq_dict[fam][gen]={}
                        if sp not in seq_dict[fam][gen]:
                            seq_dict[fam][gen][sp]={"names":[], "seqs":[]}
                        seq_dict[fam][gen][sp]["names"].append(line[1:])
                    else:
                        include=0
                else:
                    include=0
                seq=''
            else: seq +=line
        if include: seq_dict[fam][gen][sp]["seqs"].append(seq.strip().split()[0])
    return seq_dict

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def idseqdict(file):
    return_dict={}
    with open(file) as fin:
        for line in fin:
            return_dict[line.strip()]={}
    return return_dict

def filedict(file):
    return_dict={}
    with open(file) as fin:
        for line in fin:
            return_dict[line.strip()]=""
    return return_dict


def filelist(file):
    return_list=[]
    with open(file) as fin:
        for line in fin:
            return_list.append(line.strip())
    return return_list

def get_taxmap_from_file( filename ):
    return_data = {}
    if filename:
        with open( filename, 'r' ) as ranked_lineage:
            for current_line in ranked_lineage:
                key, value = process_line( current_line )
                return_data[ key ] = value
        return return_data
    return {}

def process_line( str_line ):
    split_line = str_line.split( '\t|\t' )
    return split_line[0], split_line[1]

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs
    

##########----------------------->>>>>>>>>

if __name__ == '__main__':
    main()
