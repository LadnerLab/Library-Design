#!/usr/bin/env python

# Add or remove modules here
import argparse
#import fastatools as ft		   #Available at https://github.com/jtladner/Modules
import pandas as pd

from collections import defaultdict


# the main function contains the bulk of the script
def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	# These are optional arguments
	parser.add_argument("-f", "--filler", default="........................................", help="Sequence to use as filler for designs.")
	parser.add_argument("-s", "--spacer", default=",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,", help="Sequence to use as spacer for designs.")
	
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
	
	for index,row in epi.iterrows():
			
			# Generate linear and cyclic peptides
			for item, fillLen in {"14mer": 39, "30mer": 23, "42mer": 11, "50mer": 3}.items():
				if fillLen > len(args.filler):
					print("WARNING length of user specified filler too short. Final peptides will not all be the same length.")
				epi = row[item]
				fill = args.filler[:fillLen]
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}S{epi}S", f"Linear_{item}"))
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}C{epi}C", f"Cyclic_{item}"))
			# Full length epitope, no TEV site
			epi = row["62mer"]
			pepList.append((row['Parent ID'], epi, f"S{epi}S", f"Linear_62mer"))
			pepList.append((row['Parent ID'], epi, f"C{epi}C", f"Cyclic_62mer"))
			
			# Generate the tandem peptides
			for item, spaceLen in {"14mer": 36, "30mer": 4}.items():
				if spaceLen > len(args.spacer):
					print("WARNING length of user specified spacer too short. Final peptides will not all be the same length.")
				epi = row[item]
				space = args.spacer[:spaceLen]
				pepList.append((row['Parent ID'], epi, f"{epi}{space}{epi}", f"Tandem_{item}"))
			# Generate variable length tandem peptides
			for spaceFillLen, item in {(2,25):"14mer", (16,11):"14mer", (26,1):"14mer"}.items():
				epi = row[item]
				space = args.spacer[:spaceFillLen[0]]
				fill = args.filler[:spaceFillLen[1]]
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}{epi}{space}{epi}", f"Tandem_{item}"))
			
			# Generate y-shaped peptide designs
			for spaceFillLen, item in {(0,11):"18mer", (5,1):"18mer"}.items():
				epi = row[item]
				space = args.spacer[:spaceFillLen[0]]
				fill = args.filler[:spaceFillLen[1]]
				# Regular Y-shaped
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}{epi}{space}C{THR}{epi}{space}C", f"Y-shaped_{item}"))
				# Control for Y-shaped
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}{epi}{space}S{THR}{epi}{space}S", f"Y-shaped_S-control_{item}"))
				# C-reversed Y-shaped
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}{epi}{space}C{THR}{epi[::-1]}{space}C", f"Y-shaped_C-reversed_{item}"))
				# Control for C-reversed Y-shaped
				pepList.append((row['Parent ID'], epi, f"{fill}{TEV}{epi}{space}S{THR}{epi[::-1]}{space}S", f"Y-shaped_C-reversed-S-control_{item}"))
			# Generate y-shaped peptides with no filler/TEV sites
			epi = row["18mer"]
			space = args.spacer[:10]
			pepList.append((row['Parent ID'], epi, f"{epi}{space}C{THR}{epi}{space}C", f"Y-shaped_{item}"))
			pepList.append((row['Parent ID'], epi, f"{epi}{space}S{THR}{epi}{space}S", f"Y-shaped_S-control_{item}"))
			pepList.append((row['Parent ID'], epi, f"{epi}{space}C{THR}{epi[::-1]}{space}C", f"Y-shaped_C-reversed_{item}"))
			pepList.append((row['Parent ID'], epi, f"{epi}{space}S{THR}{epi[::-1]}{space}S", f"Y-shaped_C-reversed-S-control_{item}"))

	write_tsv(pepList,args.output)

#----------------------End of main()

# Put any other function definitions here

#write output tsv file
def write_tsv(pepInfo, new_filename):
	fout=open(new_filename, 'w')
	fout.write("%s\t%s\t%s\t%s\t%s\n" % ("CodeName","Parent ID","Epitope","Peptide","Category"))
	for i in range(len(pepInfo)):
		codeName = "CT64_" + str(str(i+1).zfill(6))
		fout.write("%s\t%s\t%s\t%s\t%s\n" % (codeName,pepInfo[i][0],pepInfo[i][1],pepInfo[i][2],pepInfo[i][3]))
	fout.close()

###------------------------------------->>>>	

# This allows the script to be imported as a module to other scripts, without running the main() code
# Therefore, you could use the other defined functions in another script, if you wanted
if __name__ == "__main__":
	main()

