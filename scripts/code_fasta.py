#!/usr/bin/env python

from __future__ import division
import string, random, optparse

#This script reads in a fasta file and creates two outputs
    # 1) A new version of the fasta with names replaces by codes
    # 2) A text file linking the original names to the new coded names

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-f', '--fasta',  help='Fasta file. [None, REQ]')
    p.add_option('-l', '--length',  help='Length of code to use. If not provided, an optimal length will be calculated based on # of seqs [None]')

    opts, args = p.parse_args()
    
    #Read in fasta file
    names, seqs = read_fasta_lists(opts.fasta)
    
    #Characters to use in code
    chars=string.ascii_letters + string.digits
    
    #Determine the size of code to use, if not provided
    if not opts.length:
        length=1
        while len(chars)**length<(2*len(names)):
            length+=1
    
    #Generate codes
    codes = set()
    while len(codes)<len(names):
        codes.add(rand_code(chars,length))
#        print len(codes)
    codes=list(codes)
    
    write_fasta(codes, seqs, "coded_%s" % opts.fasta)
    write_code(names, codes, "coded_%s_key.txt" % (".".join(opts.fasta.split(".")[:-1])))
#----------------------End of main()

def write_code(names, codes, outfile):
    fout=open(outfile, "w")
    for i in range(len(names)):
        fout.write("%s\t%s\n" % (names[i], codes[i]))
    fout.close()

def rand_code(chars, l):
    code=""
    for each in range(l):
        code+=random.choice(chars)
#    print code
    return code

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

