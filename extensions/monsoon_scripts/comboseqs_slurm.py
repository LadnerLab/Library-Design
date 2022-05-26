#!/usr/bin/env python3

import argparse
from collections import defaultdict
import os
from subprocess import Popen, PIPE


parser = argparse.ArgumentParser(description='''A wrapper script for batch processing pairs of files listed in a metafile with comboseqs.py on slurm.''')

parser.add_argument("-k", "--ksize", default=7, type=int, metavar='\b', help="kmer size to use for comparison [default: 7]")
parser.add_argument("-p", "--path", default="comboseqs.py", metavar='\b', help="Path of comboseqs.py [default: comboseqs.py]")
#Resource utilization
parser.add_argument("-c", "--cpus", default=1, type=int, metavar='\b', help="Use this to specify the number of cpus per task [default: 1].")
parser.add_argument("-m", "--mem", default=1, type=int, metavar='\b', help="Use this to specify the amount of memory to allow be used, in G [default: 1].")
parser.add_argument("-t", "--time", default="20:00", metavar='\b', help="Use this to specify the amount of time needed for script to run [default: 20:00].")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-f", "--metafile", required=True, metavar='\b', help="Tab-delimited file with 2 columns, with each row containing a pair of files to be compared.")
reqArgs.add_argument("-u", "--uniprot", choices=[0,1], type=int, required=True, metavar='\b', help="Column number of UniProt files in metafile [choices: 0,1]")
reqArgs.add_argument("-g", "--genbank", choices=[0,1], type=int, required=True, metavar='\b', help="Column number of GenBank files in metafile [choices: 0,1]")

args = parser.parse_args()


#Reading in metafile and storing names of file pairs in a dictionary
metaDict = defaultdict(str)
with open(args.metafile, "r") as fin:
	for line in fin:
		line= line.strip().split("\t")
		metaDict[line[args.uniprot]]= line[args.genbank]

#Checking if directories for stdout and stderr exist; if not, create them
do= "./stdout"
de= "./stderr"
if not os.path.exists(do) and not os.path.exists(de):
	os.mkdir(do)
	os.mkdir(de)

count=0
for uniprot, genbank in metaDict.items():
	fpo= do + "/comboseqs-p" + str(count) + ".out"
	fpe= de + "/comboseqs-p" + str(count) + ".err"
	cmd= "%s -u %s -g %s -k %d" % (args.path, uniprot, genbank, args.ksize)
	cmd= "sbatch -J comboseqs-p%d -o %s -e %s -t %s -c %d --mem=%dG --wrap='%s'" % (count, fpo, fpe, args.time, args.cpus, args.mem, cmd)
	count+=1
	print (cmd)
	Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
	