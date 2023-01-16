#!/usr/bin/env python

# Add or remove modules here
import argparse
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import kmertools as kt        #Available at https://github.com/jtladner/Modules
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
    reqArgs.add_argument("-z", "--size", help="Desired size of final library.", required=True)
    reqArgs.add_argument("-o", "--output", help="Name of output TSV file.", required=True)


    args = parser.parse_args()
    
    # Add code for script here
    TEV = "GENLYFQGA"
    THR = "LVPRGS"
    
    epi = pd.read_csv(args.input,sep="\t").sort_values('Seropostive Count', ascending=False)
    epiSub = epi.iloc[:int(int(args.size) / 23)]
    #print(epi['Seropostive Count'])
    pepList = []
    category = []
    
    for index,row in epiSub.iterrows():
    	for construct in ["Linear","Cyclic","Tandem","Y-shaped","Y-shaped(S-control)","Y-shaped (C-reversed)","Y-shaped (C-reversed; S-control)"]:
    		if construct == "Linear":
    			for item in ["10mer","14mer","18mer","42mer"]:
    				category.append("Linear" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(args.filler + TEV + "S" + row[item] + "S")
    				if item == "14mer":
    					pepList.append(args.filler[:19] + TEV + "S" + row[item] + "S")
    				if item == "18mer":
    					pepList.append(args.filler[:15] + TEV + "S" + row[item] + "S")
    				if item == "42mer":
    					pepList.append("S" + row[item] + "S")
    		if construct == "Cyclic":
    			for item in ["10mer","14mer","18mer","42mer"]:
    				category.append("Cyclic" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(args.filler + TEV + "C" + row[item] + "C")
    				if item == "14mer":
    					pepList.append(args.filler[:19] + TEV + "C" + row[item] + "C")
    				if item == "18mer":
    					pepList.append(args.filler[:15] + TEV + "C" + row[item] + "C")
    				if item == "42mer":
    					pepList.append("C" + row[item] + "C")
    		if construct == "Tandem":
    			for item in ["10mer","14mer","18mer"]:
    				category.append("Tandem" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(row[item] + args.spacer + row[item])
    				if item == "14mer":
    					pepList.append(row[item] + args.spacer[:16] + row[item])
    				if item == "18mer":
    					pepList.append(row[item] + args.spacer[:8] + row[item])
    		if construct == "Y-shaped":
    			for item in ["10mer","14mer","18mer"]:
    				category.append("Y-shaped" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(row[item] + args.spacer[:8] + "C" + THR + row[item] + args.spacer[:8] + "C")
    				if item == "14mer":
    					pepList.append(row[item] + args.spacer[:4] + "C" + THR + row[item] + args.spacer[:4] + "C")
    				if item == "18mer":
    					pepList.append(row[item] + "C" + THR + row[item] + "C")
    		if construct == "Y-shaped (S-control)":
    			for item in ["10mer","14mer","18mer"]:
    				category.append("Y-shaped(S-control)" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(row[item] + args.spacer[:8] + "S" + THR + row[item] + args.spacer[:8] + "S")
    				if item == "14mer":
    					pepList.append(row[item] + args.spacer[:4] + "S" + THR + row[item] + args.spacer[:4] + "S")
    				if item == "18mer":
    					pepList.append(row[item] + "S" + THR + row[item] + "S")
    		if construct == "Y-shaped (C-reversed)":
    			for item in ["10mer","14mer","18mer"]:
    				category.append("Y-shaped (C-reversed)" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(row[item] + args.spacer[:8] + "C" + THR + row[item][::-1] + args.spacer[:8] + "C")
    				if item == "14mer":
    					pepList.append(row[item] + args.spacer[:4] + "C" + THR + row[item][::-1] + args.spacer[:4] + "C")
    				if item == "18mer":
    					pepList.append(row[item] + "C" + THR + row[item] + "C")
    		if construct == "Y-shaped (C-reversed; S-control)":
    			for item in ["10mer","14mer","18mer"]:
    				category.append("Y-shaped (C-reversed; S-control)" + "(" + item + ")")
    				if item == "10mer":
    					pepList.append(row[item] + args.spacer[:8] + "S" + THR + row[item][::-1] + args.spacer[:8] + "S")
    				if item == "14mer":
    					pepList.append(row[item] + args.spacer[:4] + "S" + THR + row[item][::-1] + args.spacer[:4] + "S")
    				if item == "18mer":
    					pepList.append(row[item] + "S" + THR + row[item][::-1] + "S")
    
    write_tsv(pepList,args.output,category)
    
#----------------------End of main()

# Put any other function definitions here

#write output tsv file
def write_tsv(peptides, new_filename,cat):
	fout=open(new_filename, 'w')
	fout.write("%s\t%s\t%s\n" % ("CodeName","Peptide","Category"))
	for i in range(len(peptides)):
		codeName = "CT44_" + str(str(i+1).zfill(6))
		fout.write("%s\t%s\t%s\n" % (codeName,peptides[i],cat[i]))
	fout.close()

###------------------------------------->>>>    

# This allows the script to be imported as a module to other scripts, without running the main() code
# Therefore, you could use the other defined functions in another script, if you wanted
if __name__ == "__main__":
    main()

