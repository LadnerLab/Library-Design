//
//  main.cpp
//  codon_sampling
//
//  Created by Paul Altin on 10.05.18.
//  Copyright © 2018 Paul Altin. All rights reserved.
//

#include <string>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <ios>
#include <iomanip>
#include <sstream>
#include <limits.h>

#include "table.h"
#include "xoroshiro.h"

#define DEBUG_INPUT 0

const uint64_t DEFAULT_TRIALS = 10000;
const uint8_t MAX_LINE_LENGTH = 128;
const uint16_t DEFAULT_NUM_SUBSAMPLE = 10000;

class FileInput
{
 public:
    std::string name;
    std::string data;

    FileInput( const char *, const char* );
    FileInput();
};

FileInput::FileInput( const char *new_name, const char *new_data )
{
    name = new_name;
    data = new_data;

}

FileInput::FileInput()
{
    name = "";
    data = "";

}


class Encoding
{
  public:
    FileInput original;

    std::string encoding;

    double gc_ratio;
    uint8_t nucleotides[ 4 ] = { 0 };
    uint8_t codons[ 64 ]     = { 0 };

    double calc_gc_ratio( void )
    {
        uint8_t g_and_c = nucleotides[ G_INDEX ] + nucleotides[ C_INDEX ];
        uint8_t t_and_a = nucleotides[ T_INDEX ] + nucleotides[ A_INDEX ];

        gc_ratio = (double) g_and_c / (double) t_and_a;

        return gc_ratio;
    }

      ~Encoding()
      {
          delete &encoding;
          delete &original;
      }
};

int encoding_compar( const void *first, const void *second )
{
    Encoding *first_ptr  = (Encoding*) first;
    Encoding *second_ptr = (Encoding*) second;

    double diff = first_ptr->gc_ratio - second_ptr->gc_ratio;

    if( diff > 0 )
        {
            return 1;
        }
    if( diff < 0 )
        {
            return -1;
        }
    return 0;
}
// custom assertion
#undef assert
#define assert(cond, message...) {  \
    if (!(cond)) {                  \
        printf(message);            \
        exit(EXIT_FAILURE);         \
    }                               \
}


const char *ARGS = "i:n:g:s:r:p:t:b:h?";
uint32_t count_lines_in_file( const char *filename );

int main(int argc, char * const argv[])
{
    // default values
    uint64_t trials = DEFAULT_TRIALS;
    uint8_t buffer = 0;
    uint16_t num_to_subsample = DEFAULT_NUM_SUBSAMPLE;
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
    uint32_t lines = count_lines_in_file( input_file );
    uint32_t loop_index = 0;
    uint32_t index = 0;
    uint32_t inner_index = 0;

    char line[MAX_LINE_LENGTH];
    std::vector<FileInput> file_data_arr;
    std::vector<std::string> results;

    file_data_arr.reserve( lines );
    results.reserve( lines );

    for( loop_index = 0; loop_index < lines; ++loop_index )
    {
        fgets( line, MAX_LINE_LENGTH, fin );
        // replace trailing newline with terminator so we can use strlen
        line[strcspn(line, "\n")] = '\0';

        // line consists of name,input
        char* input = nullptr;
        char* name = strtok_r(line, ",", &input);

        FileInput new_input( name, input );
        file_data_arr.push_back( new_input );
    }

    Encoding **encodings = (Encoding**) malloc( sizeof( Encoding *) * lines * trials );

    for( loop_index = 0; loop_index < lines; ++loop_index )
        {

            FileInput file_data = file_data_arr.at( loop_index );

            uint8_t aa_counts[ 20 ] = { 0 };
            double aa_total = file_data.data.length();
            uint64_t result_len = file_data.name.length() + 1 + digits + 1 + ( 4 * file_data.data.length() );


            for( index = 0; index < file_data.data.length(); ++index )
                {
                    ++aa_counts[ (uint8_t) acid_map[ (uint8_t) file_data.data[ index ] ] ];
                }


            uint64_t current_aa    = 0;
            uint64_t current_trial = 0;
            uint16_t len = 0;
            uint64_t i = 0;

            len = file_data.data.length();

            // trials
            for ( current_trial = 0; current_trial < trials; ++current_trial)
                {
                    Encoding *current = new Encoding();
                    // keep track of nucleotide and codon ratios
                    current->original = file_data;
            
                    // calculate result string
                    for ( current_aa = 0; current_aa < len; ++current_aa )
                        {
                            double r = xoroshiro::uniform();
                            const char aa = file_data.data[current_aa];
                           
                            codon** cod = t[aa];
                            double accum = (*cod)->w;
                           
                            while ( accum < r )
                                {
                                    accum += (*++cod)->w;
                                }
                           
                            for ( i = 0; i < 4; ++i )
                                {
                                    current->nucleotides[ i ] += (*cod)->nucleotides[i];
                                }
                            ++current->codons[(*cod)->index];
                            current->encoding.append( (*cod)->c );
            
                        }
                    current->calc_gc_ratio();
                    encodings[ ( trials * loop_index ) + current_trial ] = current;
                }
        }

    uint64_t num_encodings = trials * lines;
    // sort the encodings according to gc_ratio
    qsort( encodings, num_encodings, sizeof( Encoding *), encoding_compar );

    // find the index of the encoding with the smallest difference between low and high
    double smallest_ratio = INT_MAX;
    double ratio = 0;
    uint64_t smallest_ratio_index = 0;
    int64_t left_index = 0;
    uint64_t right_index = 0;

    std::vector<std::vector<Encoding*>> best_encodings;
    best_encodings.reserve( lines );

    for( index = 0; index < lines; index++ )
        {

            smallest_ratio = INT_MAX;
            ratio = 0;
            smallest_ratio_index = 0;
            std::vector<Encoding *> current_vector;
            current_vector.reserve( num_to_subsample );

            for( inner_index = trials * index; inner_index < ( trials * ( index + 1 ) ); inner_index++ )
                {
                    ratio = encodings[ inner_index ]->gc_ratio;
                    if( abs( ratio - gc_target_ratio ) < smallest_ratio )
                        {
                            smallest_ratio = ratio;
                            smallest_ratio_index = inner_index;
                        }
                }

            left_index  = std::max( 0LU, smallest_ratio_index - 1 );
            right_index = std::min( smallest_ratio_index + 1, num_encodings );

            current_vector.push_back( encodings[ smallest_ratio_index ] );

            while( current_vector.size() < num_to_subsample &&
                   ( left_index > 0 || right_index <= num_encodings )
                 )
                {
                    if( left_index >= 0 )
                        {
                            current_vector.push_back( encodings[ left_index ] );
                            left_index--;
                        }
                    if( right_index < num_encodings )
                        {
                            current_vector.push_back( encodings[ right_index ] );
                            right_index++;
                        }
                }
            best_encodings.push_back( current_vector );
        }

    std::string out_string;
    for( index = 0; index < lines; index++ )
        {
            std::vector<Encoding *> current_vector = best_encodings[ index ];

            for( inner_index = 0; inner_index < num_to_subsample; inner_index++ )
                {
                    std::ostringstream digit_str;

                    digit_str << std::internal << std::setfill( '0' ) << std::setw( digits + 1 ) << inner_index;
                    out_string.append( current_vector[ index ]->original.name );
                    out_string.append( "_" );

                    out_string.append( digit_str.str() );
                    out_string.append( "," );
                    out_string.append( current_vector[ index ]->original.data );
                    out_string.append( "," );
                    out_string.append( current_vector[ index ]->encoding );
                }
        }
     

    fclose(foutr);
    fclose(fouts);
    fclose(fin);

    printf("Processed %d lines x %lu trials, time elapsed %f s\n", lines, trials, omp_get_wtime() - begin );

    return EXIT_SUCCESS;
}

uint32_t count_lines_in_file( const char *filename )
{
    FILE *open_file = NULL;
    int c = 0;
    uint32_t count = 0;

    open_file = fopen( filename, "r" );

    while( ( c = fgetc( open_file ) ) != EOF )
        {
            if( !( c - (int) '\n' ) )
                {
                    count++;
                }
        }

    return count;
}
