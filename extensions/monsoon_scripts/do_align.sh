#!/bin/bash
#SBATCH --time=15:00
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH --array=0-10064
#SBATCH -o aligned.txt
#SBATCH --open-mode=append

module load muscle
module load python/3.latest

in_dir="$1"
cd "$in_dir"

files=($(ls | grep -v "out" ))
file=${files[$SLURM_ARRAY_TASK_ID]}


srun muscle -in "$file" -out "$file".out.aligned
srun protein_oligo_main.py -a "$file".out.aligned -x 9 -w 24 -s 15 -o "$file".aligned_out 

#Usage: protein_oligo_main.py [options]
#
#Options:
#  -h, --help            show this help message and exit
#  -a ALIGNMENT, --alignment=ALIGNMENT
#                        Fasta query file of sequence alignment to be used by
#                        program. [None, Required]
#  -w WINDOWSIZE, --windowSize=WINDOWSIZE
#                        Amount of characters from each alignment sequence to
#                        look at. [100]
#  -o OUTPUT, --outPut=OUTPUT
#                        Name of file program output will be written to.
#                        [oligo_out.fasta]
#  -p PERCENTVALID, --percentValid=PERCENTVALID
#                        Percent of non '-' characters present in order for the
#                        sequence to be considered valid, sequences with less
#                        than specified amount will not be present in program
#                        out put. [90.00]
#  -l MINLENGTH, --minLength=MINLENGTH
#                        Minimum length of concurrent non-dash characters that
#                        must be present in order for the sequence to be
#                        considered valid, sequences with a maximum length of
#                        concurrent non-dash characters less than this
#                        parameter will not be included in program output.
#                        [None, Required]
#  -s STEPSIZE, --stepSize=STEPSIZE
#                        Step size to move over after each subset of windowSize
#                        characters has been read
#  -x XMERWINDOWSIZE, --XmerWindowSize=XMERWINDOWSIZE
#                        Window size of Xmer sequences used in redundancy
#                        calculations [8].
#  --dont_span_gaps      Include if you do not want sliding window approach to
#                        span gaps as it walks across each sequence. Kmers with
#                        more than percentValid percent of gaps or less than
#                        minLength number of gaps will be included. Otherwise,
#                        kmer searching will be done across gaps.
