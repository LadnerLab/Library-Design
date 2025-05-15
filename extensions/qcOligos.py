#! /usr/bin/env python3

import inout as io         # Available here: https://github.com/jtladner/Modules
import fastatools as ft    # Available here: https://github.com/jtladner/Modules
import glob, optparse
from collections import defaultdict

try:
    from Bio.Seq import Seq
    testTranslation = True
except:
    testTranslation = False
    print("***Not testing translations because couldn't import biopython\n")

try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['pdf.fonttype'] = 42
    plt.rcParams["font.family"] = "Arial"
    makePlot = True
except:
    makePlot = False
    print("***Not making plot because couldn't import matplotlib\n")
    


def main():

    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
    p.add_option('-p', '--pepFasta', help='Filepath for fasta file containing designed peptides [None, Required]')
    p.add_option('-n', '--nucFasta', help='Filepath for fasta file containing encodings [None, Required]')
    p.add_option('-a', '--nucFastaAdapt', help='Filepath for fasta file containing encodings with 19mer adapter sequences added [None, Required]')
    p.add_option('-l', '--logFile', default="qcOligos.log", help='Name for log file of results [qcOligos.log]')
    p.add_option('-o', '--figName', default="gcBoxPlot.png", help='Name for log file of results [gcBoxPlot.png]')
    p.add_option('-t', '--trimTo', default=40, type=int, help='With this option, you can specify a length to trim the nucleotide encodings to before checking for uniqueness, etc. This is useful, for example, if you will not be able to seequence through the full DNA tag variable region [40]')
    
    opts, args = p.parse_args()
    
    #Open log file
    fout = open(opts.logFile, "w")
    
    if opts.pepFasta:
        fout.write("\n******Checking peptide fasta file.******\n")
        checkFasta(opts.pepFasta, fout)
        identNames(opts.pepFasta, fout)
        
        if opts.nucFasta and testTranslation:
            checkTranslation(opts.pepFasta, opts.nucFasta, fout)
    
    if opts.nucFasta:
        fout.write("\n\n******Checking encodings fasta file.******\n")
        checkFasta(opts.nucFasta, fout)
        identNames(opts.nucFasta, fout)
        findIdentical(opts.nucFasta, fout)

        if opts.trimTo:
            fout.write("\n\n******Checking encodings fasta file, trimmed to %d nucleotides.******\n" % (opts.trimTo))
            checkFasta(opts.nucFasta, fout, trimTo=opts.trimTo)
            findIdentical(opts.nucFasta, fout, trimTo=opts.trimTo)

    if opts.nucFastaAdapt:
        fout.write("\n\n******Checking encodings w/Adapters fasta file.******\n")
        checkFasta(opts.nucFastaAdapt, fout)
        checkAdapters(opts.nucFastaAdapt, fout)
        identNames(opts.nucFastaAdapt, fout)
        
        if opts.nucFasta:
            checkNucsMatch(opts.nucFasta, opts.nucFastaAdapt, fout)
    

    if makePlot:
        
        seqsD = {}
        if opts.nucFasta:
            n, s = ft.read_fasta_lists(opts.nucFasta)
            seqsD["Ecodings"] = s
        if opts.nucFastaAdapt:
            n,s = ft.read_fasta_lists(opts.nucFastaAdapt)
            seqsD["wAdapters"] = s

        if len(seqsD)>0:
            gcD = {k:[gc(x) for x in v] for k,v in seqsD.items()}
        
            fig,ax = plt.subplots(1,1,figsize=(8, 5),facecolor='w')
    
            ax.boxplot(list(gcD.values()))
            ax.set_xticks(range(1, len(gcD)+1))
            ax.set_xticklabels(gcD.keys())
            ax.set_ylabel("Proportion GC", fontsize=25)
            
            #Save figure
            fig.savefig(opts.figName,dpi=200,bbox_inches='tight')

    #Close log file
    fout.close()
    
###------------------------End of main()--------------------------------

def checkNucsMatch(ntF, adF, fout):
    ntD = ft.read_fasta_dict_upper(ntF)
    adD = ft.read_fasta_dict_upper(adF)
    
    matchCount=0
    for name, ntS in ntD.items():
        if adD[name][19:-19] != ntS:
            fout.write("%s does not match %s\n" % (ntS, adD[name][19:-19]))
        else:
            matchCount+=1
    
    if matchCount == len(ntD):
        fout.write("All encodings with adapters (%d) match the encodings without adapters.\n\n" % (matchCount))

def checkFasta(fasta, fout, trimTo=False):
    names, seqs = ft.read_fasta_lists(fasta)
    
    if trimTo:
        seqs = [s[:trimTo] for s in seqs]
    
    lengths = set([len(s) for s in seqs])
    if len(lengths)>1: 
        fout.write("There are multiple seqeunce lengths in %s: %s\n" % (fasta, ",".join([str(x) for x in lengths])))
    else:
        fout.write("All seqs in %s are %d characters long\n" % (fasta, list(lengths)[0]))
        
    uniqN = set(names)
    uniqS = set(seqs)
    
    if len(names) != len(uniqN):
        fout.write("Total Names: %d, Unique Names: %d\n" % (len(names), len(uniqN)))
    else:
        fout.write("All %d Names are Unique.\n" % (len(names)))
    if len(seqs) != len(uniqS):
        fout.write("Total Seqs: %d, Unique Seqs: %d\n\n" % (len(seqs), len(uniqS)))
    else:
        fout.write("All %d Seqs are Unique.\n\n" % (len(seqs)))


def findIdentical(fasta, fout, trimTo=False):
    names, seqs = ft.read_fasta_lists(fasta)

    if trimTo:
        seqs = [s[:trimTo] for s in seqs]

    seqD = defaultdict(list)
    
    for i,s in enumerate(seqs):
        seqD[s].append(names[i].split("-")[0])
    
    multiPeps = 0
    uniqueEnc = 0
    uniquePeps = 0
    
    for k, v in seqD.items():
        if len(v) != 1:
            if len(set(v)) != 1:
                fout.write("Non-unique seqs!  %s: %s\n" % (k, ", ".join(v)))
                multiPeps+=1
            else:
                uniquePeps+=1
        else:
            uniqueEnc+=1
    
    fout.write("%d oligos are linked to multiple peptides\n" % (multiPeps))
    fout.write("%d oligos are unique\n" % (uniqueEnc))
    fout.write("%d oligos are present multiple times for the same peptide\n" % (uniquePeps))
            
    
def identNames(fasta, fout):
    names, seqs = ft.read_fasta_lists(fasta)
    fD=defaultdict(list)
    for i in range(len(names)):
        fD[names[i]].append(seqs[i])
    for k,v in fD.items():
        if len(v)>1:
            fout.write("Mutliple identical names: %s, %d ocurrences" % (k, len(v)))
            
def checkAdapters(probes, fout):
    pN, pS = ft.read_fasta_lists(probes)
    ad = {
        "F_Adapter": "CCTATACTTCCAAGGCGCA",
        "R_Adapter": "GGTGACTCTCTGTCTTGGC",
    }
    
    fC=0
    rC=0
    
    for i in range(len(pN)):
        if not pS[i].startswith(ad["F_Adapter"]):
            fout.write("%s doesn't start with Forward adapter (%s): %s\n" % (pN[i], ad["F_Adapter"], pS[i]))
        else:
            fC+=1
            
        if not pS[i].endswith(ad["R_Adapter"]):
            fout.write("%s doesn't end with Reverse adapter (%s): %s\n" % (pN[i], ad["R_Adapter"], pS[i]))
        else:
            rC+=1
    
    if fC == len(pN):
        fout.write("All encodings (%d) have the proper forward adapter.\n\n" % (fC))
    if rC == len(pN):
        fout.write("All encodings (%d) have the proper reverse adapter.\n\n" % (rC))
    
def checkTranslation(aaF, ntF, fout):
    aaD = ft.read_fasta_dict_upper(aaF)
    ntD = ft.read_fasta_dict_upper(ntF)
    
    matchCount=0
    for name, ntS in ntD.items():
        expAA = str(Seq(ntS).translate())
        aaName = name.split("-")[0]
        if expAA != aaD[aaName]:
            fout.write("%s does not translate to %s\n" % (ntS, aaD[aaName]))
        else:
            matchCount+=1
    
    if matchCount == len(ntD):
        fout.write("All encodings (%d) match expected AA sequence.\n\n" % (matchCount))

def gc(seq):
    return (seq.count("C")+seq.count("G"))/len(seq)


###---------------------------->>>

if __name__ == "__main__":
    main()  
