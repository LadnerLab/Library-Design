//
//  main.cpp
//  codon_sampling
//
//  Created by Paul Altin on 10.05.18.
//  Copyright © 2018 Paul Altin. All rights reserved.
//

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>

#include "table.h"
#include "xoroshiro.h"

#define DEBUG_INPUT 0


// custom assertion
#undef assert
#define assert(cond, message...) {  \
    if (!(cond)) {                  \
        printf(message);            \
        exit(EXIT_FAILURE);         \
    }                               \
}


const auto MAX_LINE_LENGTH = 128;
const char *ARGS = "i:n:g:s:r:p:t:b:h?";


int main(int argc, char * const argv[])
{
    // default values
    int trials = 10000;
    uint8_t buffer = 0;
    uint16_t num_to_subsample = 0;
    double gc_target_ratio = 0;
    bool custom_buffer = false;
    const char* input_file = nullptr;
    const char* seq_output_file = nullptr;
    const char* ratio_output_file = nullptr;
    const char* probability_file = nullptr;


#if DEBUG_INPUT
    
    // debugging value
    trials = 100;
    buffer = 20;
    custom_buffer = true;
    input_file = "/Users/Altin/Desktop/artefacts/input.txt";
    seq_output_file = "/Users/Altin/Desktop/artefacts/output_sequences.txt";
    ratio_output_file = "/Users/Altin/Desktop/artefacts/output_ratios.txt";
    probability_file = "/Users/Altin/Desktop/artefacts/coding_probability_indexed.csv";

#else
    
    // parse command line options
    int opt = 0;
    while ((opt = getopt(argc, argv, ARGS )) != -1)
    {
        switch (opt)
        {
            case 'i': input_file = optarg; break;
            case 's': seq_output_file = optarg; break;
            case 'r': ratio_output_file = optarg; break;
            case 'p': probability_file = optarg; break;
            case 'n': num_to_subsample = atoi( optarg ); break;
            case 't': trials = atoi(optarg); break;
            case 'g': gc_target_ratio = atof( optarg ); break;
            case 'b': buffer = atoi(optarg); custom_buffer = true; break;
            case 'h':
            case '?':
                printf("usage: codon_sampling -i input_file -s seq_output_file -r ratio_output_file -p probability_file -n num_to_subsample -g gc_target_ratio [-t num_trials -b buffer_size]\n");
                printf("   input_file: lines must be formatted as {identifier},{sequence}, with at most %d characters per line.\n", MAX_LINE_LENGTH);
                printf("   seq_output_file: path to sequence output file (will be overwritten if it exists).\n");
                printf("   ratio_output_file: path to ratios output file (will be overwritten if it exists).\n");
                printf("   probability_file: lines must be formatted as {letter},{nucleotides,3},{weighting},{index}. The weightings do not need to sum to 1. Codon indices must range from 0 to 63.\n");
                printf("   num_to_subsample: number of 'top' encodings to take, the best encodings are measured by absolute difference from gc_target_ratio.\n" );
                printf("   gc_target_ratio: ratio GC to AT to target for encodings.\n" );
                printf("   trials: number of nucleotide sequences to generate for each input sequence. The default is 10,000.\n");
                printf("   buffer: size in bytes of output stream buffer expressed as a power of 2, between 4 and 30. Leave unspecified to use system default.\n");
                printf("Note that the optimal buffer size may depend on hardware and on the number of trials.\n");
                exit(EXIT_SUCCESS);
            default: break;
        }
    }

    // check options
    assert(input_file && probability_file, "Parameter error: need two input files. Type 'codon_sampling -h' for help.\n");
    assert(seq_output_file && ratio_output_file, "Parameter error: output file missing. Type 'codon_sampling -h' for help.\n");
    assert(trials >= 0, "Parameter error: number of trials must be non-negative. Type 'codon_sampling -h' for help.\n");
    assert(!custom_buffer || (buffer >= 4 && buffer <= 30), "Parameter error: buffer size must be between 4 and 30. Type 'codon_sampling -h' for help.\n");

#endif
    

    // start
    double begin  = omp_get_wtime();
    
    // seed RNG using hardware source
    xoroshiro::seedrandom();

    // number of digits required in suffix
    size_t digits = (size_t)floor(log10(trials)) + 1;

    // open files
    FILE* fin = fopen(input_file, "r");
    FILE* fouts = fopen(seq_output_file, "w");
    FILE* foutr = fopen(ratio_output_file, "w");
    auto t = table(probability_file);
    
    // mapping of amino acid string to index
    char acid_map['Z'] = {0};
    acid_map['A'] = 0;  acid_map['C'] = 1;  acid_map['D'] = 2;  acid_map['E'] = 3;
    acid_map['F'] = 4;  acid_map['G'] = 5;  acid_map['H'] = 6;  acid_map['I'] = 7;
    acid_map['K'] = 8;  acid_map['L'] = 9;  acid_map['M'] = 10; acid_map['N'] = 11;
    acid_map['P'] = 12; acid_map['Q'] = 13; acid_map['R'] = 14; acid_map['S'] = 15;
    acid_map['T'] = 16; acid_map['V'] = 17; acid_map['W'] = 18; acid_map['Y'] = 19;

    // set size of stream buffer if requested (use full buffering)
    if (custom_buffer)
        {
            setvbuf(fouts, nullptr, _IOFBF, 1 << buffer);
            setvbuf(foutr, nullptr, _IOFBF, 1 << buffer);
        }
    
    // process line by line
    int lines = 0;
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fin))
    {
        ++lines;

        // replace trailing newline with terminator so we can use strlen
        line[strcspn(line, "\n")] = 0;

        // line consists of name,input
        char* input = nullptr;
        char* name = strtok_r(line, ",", &input);
        
        // length of each component
        size_t len = strlen(input);
        size_t namelen = strlen(name);
        size_t inputlen = strlen(input);
        size_t resultlen = namelen+1+digits+1 + inputlen+1 + 3*len;
        
        uint16_t i = 0;
        uint16_t j = 0;

        // calculate amino acid ratios
        size_t aa_counts[20] = {0};
        for ( i = 0; i < len; ++i)
            ++aa_counts[ (uint8_t) acid_map[ (uint8_t) input[i] ] ];
        double aa_total = len;
        
        // prepare an array to store result (fwrite is much faster than fprintf)
        char result[resultlen+1];
        result[resultlen] = '\n';

        memcpy(&result[0], name, namelen);
        result[namelen] = '_';
        for ( i = 0; i < digits; ++i)
            {
                result[namelen+1+i] = '0';
            }

        result[namelen+1+digits] = ',';
        memcpy(&result[namelen+1+digits+1], input, inputlen);
        result[namelen+1+digits+1+inputlen] = ',';
        
        // trials
        for ( j = 0; j < trials; ++j)
        {
            // keep track of nucleotide and codon ratios
            size_t nucleotides[4] = {0};
            size_t codons[64] = {0};

            // calculate result string
            for ( i = 0; i < len; ++i)
            {
                double r = xoroshiro::uniform();
                const char aa = input[i];
                
                codon** cod = t[aa];
                double accum = (*cod)->w;
                
                while (accum < r)
                    {
                        accum += (*++cod)->w;
                    }
                
                for (int j = 0; j < 4; ++j)
                    nucleotides[j] += (*cod)->nucleotides[j];
                
                ++codons[(*cod)->index];

                memcpy( &result[ namelen+1+digits+1+inputlen+1+3*i ], (*cod)->c, 3);
            }
            
            // add suffix to name (sprintf is too slow)
            for ( i = 0; i < digits; ++i)
            {
                char* d = &result[namelen+1+digits-1-i];
                if (*d == '9')
                    {
                        *d = '0';
                    }
                else
                    {
                        ++*d;
                        break;
                    }
            }
            
            // write line to file
            fwrite(result, resultlen+1, 1, fouts);
            
            // write ratios to second file
            size_t total_nucleotides = 3*len, total_codons = len;
            
            fprintf(foutr, "%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g\n",
                    (double)nucleotides[0]/total_nucleotides,
                    (double)nucleotides[1]/total_nucleotides,
                    (double)nucleotides[2]/total_nucleotides,
                    (double)nucleotides[3]/total_nucleotides,
                    (double)aa_counts[0]/aa_total,
                    (double)aa_counts[1]/aa_total,
                    (double)aa_counts[2]/aa_total,
                    (double)aa_counts[3]/aa_total,
                    (double)aa_counts[4]/aa_total,
                    (double)aa_counts[5]/aa_total,
                    (double)aa_counts[6]/aa_total,
                    (double)aa_counts[7]/aa_total,
                    (double)aa_counts[8]/aa_total,
                    (double)aa_counts[9]/aa_total,
                    (double)aa_counts[10]/aa_total,
                    (double)aa_counts[11]/aa_total,
                    (double)aa_counts[12]/aa_total,
                    (double)aa_counts[13]/aa_total,
                    (double)aa_counts[14]/aa_total,
                    (double)aa_counts[15]/aa_total,
                    (double)aa_counts[16]/aa_total,
                    (double)aa_counts[17]/aa_total,
                    (double)aa_counts[18]/aa_total,
                    (double)aa_counts[19]/aa_total,
                    (double)codons[0]/total_codons,
                    (double)codons[1]/total_codons,
                    (double)codons[2]/total_codons,
                    (double)codons[3]/total_codons,
                    (double)codons[4]/total_codons,
                    (double)codons[5]/total_codons,
                    (double)codons[6]/total_codons,
                    (double)codons[7]/total_codons,
                    (double)codons[8]/total_codons,
                    (double)codons[9]/total_codons,
                    (double)codons[10]/total_codons,
                    (double)codons[11]/total_codons,
                    (double)codons[12]/total_codons,
                    (double)codons[13]/total_codons,
                    (double)codons[14]/total_codons,
                    (double)codons[15]/total_codons,
                    (double)codons[16]/total_codons,
                    (double)codons[17]/total_codons,
                    (double)codons[18]/total_codons,
                    (double)codons[19]/total_codons,
                    (double)codons[20]/total_codons,
                    (double)codons[21]/total_codons,
                    (double)codons[22]/total_codons,
                    (double)codons[23]/total_codons,
                    (double)codons[24]/total_codons,
                    (double)codons[25]/total_codons,
                    (double)codons[26]/total_codons,
                    (double)codons[27]/total_codons,
                    (double)codons[28]/total_codons,
                    (double)codons[29]/total_codons,
                    (double)codons[30]/total_codons,
                    (double)codons[31]/total_codons,
                    (double)codons[32]/total_codons,
                    (double)codons[33]/total_codons,
                    (double)codons[34]/total_codons,
                    (double)codons[35]/total_codons,
                    (double)codons[36]/total_codons,
                    (double)codons[37]/total_codons,
                    (double)codons[38]/total_codons,
                    (double)codons[39]/total_codons,
                    (double)codons[40]/total_codons,
                    (double)codons[41]/total_codons,
                    (double)codons[42]/total_codons,
                    (double)codons[43]/total_codons,
                    (double)codons[44]/total_codons,
                    (double)codons[45]/total_codons,
                    (double)codons[46]/total_codons,
                    (double)codons[47]/total_codons,
                    (double)codons[48]/total_codons,
                    (double)codons[49]/total_codons,
                    (double)codons[50]/total_codons,
                    (double)codons[51]/total_codons,
                    (double)codons[52]/total_codons,
                    (double)codons[53]/total_codons,
                    (double)codons[54]/total_codons,
                    (double)codons[55]/total_codons,
                    (double)codons[56]/total_codons,
                    (double)codons[57]/total_codons,
                    (double)codons[58]/total_codons,
                    (double)codons[59]/total_codons,
                    (double)codons[60]/total_codons,
                    (double)codons[61]/total_codons,
                    (double)codons[62]/total_codons,
                    (double)codons[63]/total_codons);
        }
    }

    fclose(foutr);
    fclose(fouts);
    fclose(fin);

    printf("Processed %d lines x %d trials, time elapsed %f s\n", lines, trials, omp_get_wtime() - begin );

    return EXIT_SUCCESS;
}
