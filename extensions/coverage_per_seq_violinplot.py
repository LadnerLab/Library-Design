#!/usr/bin/env python3

import argparse, re
import kmertools as kt		#Available at https://github.com/jtladner/Modules
import fastatools as ft		#Available at https://github.com/jtladner/Modules
import pandas as pd
import seaborn as sns

#Example command: coverage_per_seq_violinplot.py -d /Users/colleenung/Documents/197911_InfluenzavirusA/HA/SW_SC_noC/t0.200/197911_id_70_9_SWSC-x9-y30-t0.200.fasta -c /Users/colleenung/Documents/197911_InfluenzavirusA/HA/197911_id_70_9 -k 9 -t 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95 --swCtoS -o 197911_id_70_9_coverage_per_seq_violinplot.png -s 197911_id_70_9_coverage_per_seq_stats.txt

parser = argparse.ArgumentParser(description='''A script that will generate violin plot(s) to visualize the distribution of kmer coverage 
						in the design on a per sequence basis. Can generate multiple violin plots, with each representing a different Xmer 
						threshold.''')

parser.add_argument("-k", "--ksize", default=9, type=int, metavar='\b', help="Size of kmer to use for looking at kmer coverage in the design [default: 9].")
parser.add_argument("-t", "--targets", default="0.5,0.75,1", metavar='\b', help="Target thresholds to generate violin plots for.")
parser.add_argument("-o", "--output", default="coverage_per_seq_violinplot.png", metavar='\b', help="Name of output PNG file with violin plot(s). [default: coverage_per_seq_violinplot.png]")
parser.add_argument("-s", "--statsoutput", default="coverage_per_seq_violinplot.txt", metavar='\b', help="Name of output txt file with descriptive statistics. [default: coverage_per_seq_violinplot.txt]")
parser.add_argument("--swCtoS", default=False, action="store_true", help="Use this flag if Cysteine residues were converted to Serine residues in the SW portion of the design.")

#New argument group to underscore that these arguments are required despite being provided with flags
reqArgs = parser.add_argument_group("required arguments")
reqArgs.add_argument("-d", "--design", required=True, metavar='\b', help="Input design file. If looking at multiple Xmer thresholds, only provide path to one of the design files. Assuming designs share the same naming structure and are located in a directory containing subdirectories for each Xmer target threshold.")
reqArgs.add_argument("-c", "--cluster", required=True, metavar='\b', help="Input cluster file to look at kmer coverage on a per sequence basis. Note, cluster names must end with cluster number.")

args = parser.parse_args()


#Parsing target thresholds
targetThresh = sorted(list(set([float(x) for x in args.targets.split(",")])))

#Reading in fasta file (in this case, cluster file). Returns two lists, the first containing seq names and the second containing its sequences.
names, seqs = ft.read_fasta_lists(args.cluster)

xthrList=[]
coverageperseqList=[]
for thr in targetThresh:
	#Using path of input design file to find design files for other desired target threshold(s), if applicable
	searchstr= ".*/t([\d.]*)/.*"
	regexresult= re.search(searchstr, args.design)
	designPath= re.sub(str(regexresult.group(1)), ("%.3f" % (thr)), args.design)

	#Creating set of all unique kmers within design
	designkSet= kt.kmerSetFasta(designPath, args.ksize, filter=[])

	for s in seqs:
		if args.swCtoS:
			s = s.replace("C", "S")
		#Creating set of all unique kmers within sequence
		sSet = kt.kmerSet(s, args.ksize, filter=["X"])
		
		xmersCovered= sSet.intersection(designkSet)
		percentCovered= (len(xmersCovered) / len(sSet))*100
		xthrList.append(("%.3f" % (thr)))
		coverageperseqList.append(percentCovered)

labelY= "%% %dmers covered per sequence" % args.ksize
dataDict= {"Xmer Threshold":xthrList, labelY:coverageperseqList}
#Creating pandas dataframe from dictionary
df = pd.DataFrame(dataDict)


#Generating violin plot from pandas dataframe using Seaborn
ax = sns.violinplot(x=df["Xmer Threshold"], y=df[labelY], palette="Set3")
ax.set_ylabel(labelY)
ax.set_xlabel("Xmer Threshold")
fig = ax.get_figure()
fig.savefig(args.output, bbox_inches='tight', dpi=200)


#Writing out file with descriptive statistics
with open(args.statsoutput, "w") as fout:
	line1= "\tMaximum\tQ3\tMedian\tQ1\tMinimum\tIQR"
	fout.write(line1)
	
	for thr in targetThresh:
		thrDF= df.loc[df["Xmer Threshold"] == ("%.3f" % (thr))]
		
		maximum= thrDF[labelY].max()
		q3= thrDF[labelY].quantile(q=0.75, interpolation='midpoint')
		median= thrDF[labelY].quantile(q=0.5, interpolation='midpoint')
		q1= thrDF[labelY].quantile(q=0.25, interpolation='midpoint')
		minimum= thrDF[labelY].min()
		IQR= q3-q1
		
		line2= "\n%.3f\t%f\t%f\t%f\t%f\t%f\t%f" % (thr,maximum,q3,median,q1,minimum,IQR)
		fout.write(line2)