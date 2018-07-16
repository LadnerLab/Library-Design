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
    p.add_option('-k', '--kmers', default="5,6,7,8,9,10,11,12", help='Comma-delimited string of kmer sizes to check. [5,6,7,8,9,10,11,12]')
    p.add_option('--makemap', default=False, action = "store_true", help='Use this flag if you want to generate a map linking oligos to the ref seqs containing matching epitopes. Only the smallest kmers size provided is used to generate map [False]')

    opts, args = p.parse_args()
    
    #Parse kmer sizes from input
    ksizes = [int(x) for x in opts.kmers.split(",")]
    
    #Run some basic tests on designed oligos
    oligo_dict_sets = basic_tests(opts)
    
    print "kmerSize\tRef#\tDesign#\t#RefOnly\t#DesignOnly\t%RefOnly\t%DesignOnly\tAvgRedundancy\t%RefInDesign"
    
    for k in ksizes:
        ref_k_dict = fasta2kmers(opts.ref, k)
        design_k_dict = fasta2kmers(opts.oligos, k)
        
        total_ref_ks = combine_dicts(ref_k_dict.values())
        total_design_ks = combine_dicts(design_k_dict.values())
        
        ref_ks = set(total_ref_ks)
        design_ks = set(total_design_ks)
        
        in_design_only = design_ks.difference(ref_ks)
        in_ref_only = ref_ks.difference(design_ks)
        print "%d\t%d\t%d\t%d\t%d\t%.3f%%\t%.3f%%\t%.3f\t%.3f%%" % (k, len(ref_ks), len(design_ks), len(in_ref_only), len(in_design_only), len(in_ref_only)/len(ref_ks)*100, len(in_design_only)/len(design_ks)*100, len(total_design_ks)/len(design_ks), len(design_ks)/len(ref_ks)*100)
#        print "%d: %.2f%% in Ref only, %.2f%% in Design only" % (k, len(in_ref_only)/len(ref_k_dict)*100, len(in_design_only)/len(design_k_dict)*100)
        
        if opts.makemap and k == min(ksizes):
            kmer_ref_dict = {}
            for r, ks in ref_k_dict.iteritems():
                for each in ks:
                    if each not in kmer_ref_dict: kmer_ref_dict[each] = [r]
                    else: kmer_ref_dict[each].append(r)
            for d, ks in design_k_dict.iteritems():
                for each in ks:
                    oligo_dict_sets[d].update(kmer_ref_dict[each])

            fout = open("%s_epitopemap.txt" % opts.out, "w")
            for oname, rmatches in oligo_dict_sets.iteritems():
                fout.write("%s\t%s\n" % (oname, ",".join(rmatches)))

#----------------------End of main()

def basic_tests(opts):
    oligo_dict = read_fasta_dict_upper(opts.oligos)
    # Check that all oligos are the same length
    oligo_lengths = [len(x) for x in oligo_dict.values()]
    if len(set(oligo_lengths)) > 1: print "!!!!! Oligos are not all the same length. The following lengths were observed: %s" % (", ".join([str(x) for x in set(oligo_lengths)]))
    else: 
#        print "All oligos have %d amino acids" % (oligo_lengths[0])
        # Check that all oligos are present with the reference set of proteins
        ref_y_dict = fasta2kmers(opts.ref, oligo_lengths[0])

        total_ref_ys = combine_dicts(ref_y_dict.values())
        total_design_ys = oligo_dict.values()

        ref_ys = set(total_ref_ys)
        design_ys = set(total_design_ys)

        # Check that all oligos are unique
        if len(total_design_ys) != len(design_ys): print "!!!!! Out of %d oligos, only %d are unique" % (total_design_ys, design_ys)

        in_design_only = design_ys.difference(ref_ys)
        in_ref_only = ref_ys.difference(design_ys)

        if in_design_only: print "!!!!!! There are %d oligos in the design that are not in the reference" % in_design_only

        print "OligoSize\tRef#\tDesign#\t#RefOnly\t#DesignOnly\t%RefOnly\t%DesignOnly\tAvgRedundancy\t%RefInDesign"
        print "%d\t%d\t%d\t%d\t%d\t%.3f%%\t%.3f%%\t%.3f\t%.3f%%" % (oligo_lengths[0], len(ref_ys), len(design_ys), len(in_ref_only), len(in_design_only), len(in_ref_only)/len(ref_ys)*100, len(in_design_only)/len(design_ys)*100, len(total_design_ys)/len(design_ys), len(design_ys)/len(ref_ys)*100)
    return {x:set([]) for x in oligo_dict}

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
                seqs.append(seq)
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

