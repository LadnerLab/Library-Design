#!/usr/bin/env python

# Add or remove modules here
import argparse
import fastatools as ft		   #Available at https://github.com/jtladner/Modules
import kmertools as kt		  #Available at https://github.com/jtladner/Modules
import glob, os, shutil
import pandas as pd

from collections import defaultdict


# the main function contains the bulk of the script
def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	# These are optional arguments
	parser.add_argument("-f", "--filler", default="XXXXXXXXXXXXXXXXXXXXXXX", help="Sequence to use as filler for designs.")
	parser.add_argument("-s", "--spacer", default="XXXXXXXXXXXXXXXXXXXXXXXX", help="Sequence to use as spacer for designs.")
	
	# These are arguments that the user is required to provide
	reqArgs = parser.add_argument_group('required arguments')
	reqArgs.add_argument("-i", "--input", help="Base sequences for design. Tab-delimited file with ... .", required=True)
	reqArgs.add_argument("-o", "--output", help="Name of output TSV file.", required=True)


	args = parser.parse_args()
	
	# Add code for script here
	TEV = "GENLYFQGA"
	THR = "LVPRGS"
	
	epi = pd.read_csv(args.input,sep="\t").sort_values('Seropostive Count', ascending=False)
	pepList = []
	category = []
	
	for index,row in epi.iterrows():
			for item, fillLen in {"10mer": 23, "14mer": 19, "18mer": 15,"42mer": 0}.items():
				category = "Linear" + "(" + item + ")"
				if item == "42mer":
					pepList.append((row['Parent ID'],row[item],"S" + row[item] + "S",category))
				else:
					pepList.append((row['Parent ID'],row[item],args.filler[:fillLen] + TEV + "S" + row[item] + "S",category))
				print(pepList)
			for item, fillLen in {"10mer": 23, "14mer": 19, "18mer": 15,"42mer": 0}.items():
				category = "Cyclic" + "(" + item + ")"
				if item == "42mer":
					pepList.append((row['Parent ID'],row[item],"C" + row[item] + "C",category))
				else:
					pepList.append((row['Parent ID'],row[item],args.filler[:fillLen] + TEV + "C" + row[item] + "C",category))
			for item, spaceLen in {"10mer": 24, "14mer": 16, "18mer": 8}.items():
				category = "Tandem" + "(" + item + ")"
				pepList.append((row['Parent ID'],row[item],row[item] + args.spacer[:spaceLen] + row[item],category))
			for item, spaceLen in {"10mer": 8, "14mer": 4, "18mer": 0}.items():
				category = "Y-shaped" + "(" + item + ")"
				pepList.append((row['Parent ID'],row[item],row[item] + args.spacer[:spaceLen] + "C" + THR + row[item] + args.spacer[:spaceLen] + "C",category))
			for item, spaceLen in {"10mer": 8, "14mer": 4, "18mer": 0}.items():
				category = "Y-shaped(S-control)" + "(" + item + ")"
				pepList.append((row['Parent ID'],row[item],row[item] + args.spacer[:spaceLen] + "S" + THR + row[item] + args.spacer[:spaceLen] + "S",category))
			for item, spaceLen in {"10mer": 8, "14mer": 4, "18mer": 0}.items():
				category = "Y-shaped (C-reversed)" + "(" + item + ")"
				pepList.append((row['Parent ID'],row[item],row[item] + args.spacer[:spaceLen] + "C" + THR + row[item][::-1] + args.spacer[:spaceLen] + "C",category))
			for item, spaceLen in {"10mer": 8, "14mer": 4, "18mer": 0}.items():
				category = "Y-shaped (C-reversed; S-control)" + "(" + item + ")"
				pepList.append((row['Parent ID'],row[item],row[item] + args.spacer[:spaceLen] + "S" + THR + row[item][::-1] + args.spacer[:spaceLen] + "S",category))

	
	write_tsv(pepList,args.output)
	
#----------------------End of main()

# Put any other function definitions here

#write output tsv file
def write_tsv(pepInfo, new_filename):
	fout=open(new_filename, 'w')
	fout.write("%s\t%s\t%s\t%s\t%s\n" % ("CodeName","Parent ID","Epitope","Peptide","Category"))
	for i in range(len(pepInfo)):
		codeName = "CT44_" + str(str(i+1).zfill(6))
		fout.write("%s\t%s\t%s\t%s\t%s\n" % (codeName,pepInfo[i][0],pepInfo[i][1],pepInfo[i][2],pepInfo[i][3]))
	fout.close()

###------------------------------------->>>>	

# This allows the script to be imported as a module to other scripts, without running the main() code
# Therefore, you could use the other defined functions in another script, if you wanted
if __name__ == "__main__":
	main()

