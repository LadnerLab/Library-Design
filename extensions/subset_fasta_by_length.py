#! /usr/bin/env python

import sys, optparse

p = optparse.OptionParser()

p.add_option('-f', '--fastaFile', help='Input fasta file')
p.add_option('-o', '--outputName', help='Name for output fasta file')
p.add_option('-s', '--minSize', type='int', help='An integer specifying the minimum acceptable size of a sequence to be ouput')

options, args = p.parse_args()

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


###------------------>>>

if __name__ == '__main__':
	sub_names = []
	sub_seqs = []
	names, seqs = read_fasta_lists(options.fastaFile)
	for i in range(len(seqs)):
		if len(seqs[i]) >= options.minSize: 
			sub_names.append(names[i])
			sub_seqs.append(seqs[i])
	write_fasta(sub_names, sub_seqs, options.outputName)
	