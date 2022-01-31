#!/usr/bin/env python

from __future__ import division
import optparse, os
#from subprocess import Popen, PIPE
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

#This script reads in a fasta file of designed oligos and a fasta file of the protein sequences used to design the oligos
    #It then compares these sequences to verify that the designed oligos are as expected

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-g', '--oligos',  help='Fasta file with designed oligos. [None, REQ]')
    p.add_option('-r', '--ref',  help='Fasta file with reference protein sequences. [None, REQ]')
    p.add_option('-o', '--out',  help='base name for output files. [None, REQ]')
    p.add_option('-k', '--kmer', type='int', default="7", help='kmer size to use for map. [7]')
    p.add_option('--makemap', default=False, action = "store_true", help='Use this flag if you want to generate a map linking oligos to the ref seqs containing matching epitopes. Only the smallest kmers size provided is used to generate map [False]')

    opts, args = p.parse_args()
    
    #Make dict with keys=kmers and values=all reference names
    ref_k_dict = fasta2kmers(opts.ref, opts.kmer)
    kmer_ref_dict = {}
    for r, ks in ref_k_dict.iteritems():
        for each in ks:
            if each not in kmer_ref_dict: kmer_ref_dict[each] = [r]
            else: kmer_ref_dict[each].append(r)

#    design_k_dict = fasta2kmers(opts.oligos, k)
    fout = open("%s_epitopemap.txt" % opts.out, "w")

    oligo_fasta_dict = read_fasta_dict_upper(opts.oligos)
    for name, seq in oligo_fasta_dict.iteritems():
        matches=set()
        for each in seq2kmers(seq,opts.kmer):
            if each in kmer_ref_dict:
                matches.update(kmer_ref_dict[each])
            else: print "%s not in ref seqs" % each
        fout.write("%s\t%s\n" % (name, "~".join(matches)))

    fout.close()

#----------------------End of main()


#def basic_tests(opts):
#    oligo_dict = read_fasta_dict_upper(opts.oligos)
#    # Check that all oligos are the same length
#    oligo_lengths = [len(x) for x in oligo_dict.values()]
#    if len(set(oligo_lengths)) > 1: print "!!!!! Oligos are not all the same length. The following lengths were observed: %s" % (", ".join([str(x) for x in set(oligo_lengths)]))
#     else: 
# #        print "All oligos have %d amino acids" % (oligo_lengths[0])
#         # Check that all oligos are present with the reference set of proteins
#         ref_y_dict = fasta2kmers(opts.ref, oligo_lengths[0])
# 
#         total_ref_ys = combine_dicts(ref_y_dict.values())
#         total_design_ys = oligo_dict.values()
# 
#         ref_ys = set(total_ref_ys)
#         design_ys = set(total_design_ys)
# 
#         # Check that all oligos are unique
#         if len(total_design_ys) != len(design_ys): print "!!!!! Out of %d oligos, only %d are unique" % (total_design_ys, design_ys)
# 
#         in_design_only = design_ys.difference(ref_ys)
#         in_ref_only = ref_ys.difference(design_ys)
# 
#         if in_design_only: print "!!!!!! There are %d oligos in the design that are not in the reference" % in_design_only
# 
#         print "OligoSize\tRef#\tDesign#\t#RefOnly\t#DesignOnly\t%RefOnly\t%DesignOnly\tAvgRedundancy\t%RefInDesign"
#         print "%d\t%d\t%d\t%d\t%d\t%.3f%%\t%.3f%%\t%.3f\t%.3f%%" % (oligo_lengths[0], len(ref_ys), len(design_ys), len(in_ref_only), len(in_design_only), len(in_ref_only)/len(ref_ys)*100, len(in_design_only)/len(design_ys)*100, len(total_design_ys)/len(design_ys), len(design_ys)/len(ref_ys)*100)
#    return {x:set([]) for x in oligo_dict}

def combine_dicts(list_of_dicts):
    combolist = []
    for l in list_of_dicts:
        combolist+=l.keys()
    return combolist


def combine_lists(list_of_lists):
    combolist = []
    for l in list_of_lists:
        combolist+=l
    return combolist

def fasta2kmers(fasta, kmersize):
    kmer_dict = {}
    fasta_dict = read_fasta_dict_upper(fasta)
    for name, seq in fasta_dict.iteritems():
        kmer_dict[name] = seq2kmers(seq, kmersize)
    return kmer_dict

def seq2kmers(seq, kmersize):
    kmers={}
    i=0
    while i+kmersize <= len(seq):
        k = seq[i:i+kmersize]
        if not set(['-', 'X']).intersection(set([x for x in k])):
            kmers[k] = ""
        i+=1
    return kmers

def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

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
                # Only grab the first item in the line,
                # protects against SEQ followed by anything else
                seqs.append( seq.strip().split()[ 0 ] )
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()



###------------------------------------->>>>    

if __name__ == "__main__":
    main()

