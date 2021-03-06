# Nucleotide Encoding of Amino Acid Oligonucleotides
Scripts to encode a designed library of amino acid oligonupeptides as nucleotide oligonucleoties. 

### Dependencies/Requirements
    - If you plan to use the provided DeepLearning model, h2o 3.20.0.8 MUST be installed on your system.
    - OpenMP 
    - g++
    - Python 3.5 or greater
    - Pandas


### Compiling
In order to generate encodings for your protein sequences, the Cpp files must
be compiled. To do this:
```
    make optimized
```
This produces the executable called 'main', which can be called directly or can be called
by encoding_with_nn.py.

### Example 1: Creating random encodings and selecting the best ones in separate steps.
#### Step 1: 
To create 10,000 encodings for each sequence, but only output the 300 encodings that have the
lowest absolute deviation of gc content ratio from 0.55. We also use two cores for this analysis,
as denoted by the '-c 2' option. input_file must contain lines of the form {seq},{name}, where
the length of each line can be a maximum of 128. After completion, output_ratio will contain the
necessary information to input to the neural network, out_seqs will contain the top 300 encodings for
each sequence.
    
```
    ./main -i input_file -r output_ratio -p codon_weights.csv -s out_seqs -n 300 -c 2 -t 10000 -g 0.55
```
#### Step 2:

Either using the previously created 'output_reatio' and 'out_seqs' files, use deeplearning_model
to predict the best encodings in out_seqs using the data in output_ratio. In this example,
10 sequences (of k encodings each) will be processed by the Neural Network at a time. The
--read_per_loop flag should be lower of machines with less memory, and higher on machines with more.
Only the encoding with the lowest absolute neural network prediction will be output for each input sequence.
```
./encoding_with_nn.py -m deeplearning_model -r output_ratio -s out_seqs -o best_encodings --read_per_loop 10 -n 1
```
### Example 1: Creating random encodings and selecting the best ones in one step.
We can also randomly generate encodings and evaluate them with the neural network
in one step. Note that this command is equivalent to calling the two commands in step 1 in
successsion. 
```
./encoding_with_nn.py -m deeplearning_model -r output_ratio -s out_seqs -o best_encodings
    --subsample 300 --read_per_loop 10 -n 1 -c 2 -p codon_weights.csv -i input_file -t 10000
```
### Usage: main (C++)
```
 usage: codon_sampling -i input_file -s seq_output_file -r ratio_output_file -p probability_file -n num_to_subsample -g gc_target_ratio [-t num_trials ]
   input_file: lines must be formatted as {identifier},{sequence}, with at most 128 characters per line.
   seq_output_file: path to sequence output file (will be overwritten if it exists).
   ratio_output_file: path to ratios output file (will be overwritten if it exists).
   probability_file: lines must be formatted as {letter},{nucleotides,3},{weighting},{index}. The weightings do not need to sum to 1. Codon indices must range from 0 to 63.
   num_to_subsample: number of 'top' encodings to take, the best encodings are measured by absolute difference from gc_target_ratio.
   number of threads to use for operations, default is 1
   gc_target_ratio: ratio GC to AT to target for encodings.
   trials: number of nucleotide sequences to generate for each input sequence. The default is 10,000.
``` 
### Usage: encoding_with_nn.py
```
Use h2o to select encodings for oligos.

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL, --model MODEL
                        Trained model that ill be used to predict bindings
  -r RATIO_FILE, --ratio_file RATIO_FILE
                        File containing output produced by oligo_encoding
                        script.
  -n NN_SUBSET_SIZE, --nn_subset_size NN_SUBSET_SIZE
                        The nn_subset_size number of encodings will be output
                        that have the smallest nn deviation
  -o OUT_FILE, --out_file OUT_FILE
                        File to write final encodings out to.
  -s SEQUENCES, --sequences SEQUENCES
                        File to read sequences generated by the oligo_encoding
                        script
  -i INPUT, --input INPUT
                        Input file containing name,seq pairs, the maximum line
                        length is 128
  -c CORES, --cores CORES
                        Number of cores to use in the oligo_encoding script,
                        default is 2, as 2 generally has the best performance
  --subsample SUBSAMPLE
                        Number of randomly generated encodings to subsample
  --gc_target GC_TARGET
                        Target gc ratio for generated seqsThe first subsample
                        number of encodings to take will be those with minimum
                        absolute value difference from gc_target.
  -p PROBABILITY_FILE, --probability_file PROBABILITY_FILE
                        probability_file: lines must be formatted as
                        {letter},{nucleotides,3},{weighting},{index}. The
                        weightings do not need to sum to 1. Codon indices must
                        range from 0 to 63.
  -t TRIALS, --trials TRIALS
                        Number of trials to perform. For each sequence in
                        sequences, trials number of candidate encodings will
                        be created.
  --read_per_loop READ_PER_LOOP
                        Number of lines from the output seq and ratio files to
                        read at a time, the higher this parameter is the more
                        memory will be used by h2o.


```	