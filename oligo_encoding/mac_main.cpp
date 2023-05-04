
//


#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <ios>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <limits.h>
#include <unordered_set>
#include <string>
#include <thread>
#include <mutex>

#include "table.h"
#include "xoroshiro.h"

#define DEBUG_INPUT 0

const uint64_t DEFAULT_TRIALS = 10000;
const uint16_t DEFAULT_NUM_SUBSAMPLE = 10000;
const uint8_t NUM_AMINO_ACIDS = 20;
const uint8_t NUM_CODONS_POSSIBLE = 64;
const uint8_t NUM_NUCLEOTIDES = 4;
const uint8_t DEFAULT_NUM_THREADS = 1;
const uint16_t MAX_LENGTH = 128;

class FileInput {
public:
    std::string name;
    std::string data;

    uint32_t aa_total;
    uint8_t aa_counts[ 20 ] = {0};
    uint64_t total_nucleotides = 0;

    FileInput( const char *, const char * );
    FileInput();
};

FileInput::FileInput( const char *new_name, const char *new_data ) {
    name = new_name;
    data = new_data;
}

FileInput::FileInput() {
    name = "";
    data = "";
}

int valid_line( char *input_line, uint16_t max_line_length );
void report_bad_line( char *line, int lineno );

class Encoding {
  public:
    FileInput original;

    std::string encoding;

    long double gc_ratio = 0;
    long double gc_dist_abs = 0;
    uint64_t nucleotides[ 4 ] = {0};
    uint64_t codons[ 64 ]     = {0};
    uint64_t total_codons = 0;

    double calc_gc_ratio( void ) {
        uint64_t g_and_c = nucleotides[ G_INDEX ] + nucleotides[ C_INDEX ];
        uint64_t t_and_a = nucleotides[ T_INDEX ] + nucleotides[ A_INDEX ] + g_and_c;
        gc_ratio = (long double) g_and_c / (long double) t_and_a;
        return gc_ratio;
    }
};

int encoding_compar( const void *first, const void *second ) {
    Encoding **first_ptr  = (Encoding**) first;
    Encoding **second_ptr = (Encoding**) second;

    long double diff =  (*first_ptr)->gc_dist_abs - (*second_ptr)->gc_dist_abs;

    if( diff > 0 )
        {return 1;
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

const char *ARGS = "i:c:n:g:s:r:p:t:b:h:l";
uint64_t count_lines_in_file( const char *filename );

int main(int argc, char * const argv[])
{
    // default values
    uint64_t trials = DEFAULT_TRIALS;
    uint16_t num_to_subsample = DEFAULT_NUM_SUBSAMPLE;
    uint8_t num_threads = DEFAULT_NUM_THREADS;
    uint16_t max_line_length = MAX_LENGTH;
    double gc_target_ratio = 0.55;
    const char* input_file = nullptr;
    const char* seq_output_file = nullptr;
    const char* ratio_output_file = nullptr;
    const char* probability_file = nullptr;


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
                case 'c': num_threads = atoi( optarg ); break;
                case 't': trials = atoi(optarg); break;
                case 'g': gc_target_ratio = atof( optarg ); break;
                case 'h': 
                case 'l': max_line_length = atoi( optarg ) + 1; break;
                case '?':
                    printf("usage: codon_sampling -i input_file -s seq_output_file -r ratio_output_file -p probability_file -n num_to_subsample -g gc_target_ratio [-t num_trials ] -l max_line_length\n");
                    printf("   input_file: lines must be formatted as {identifier},{sequence}, with at most 65,534 characters per line.\n");
                    printf("   seq_output_file: path to sequence output file (will be overwritten if it exists).\n");
                    printf("   ratio_output_file: path to ratios output file (will be overwritten if it exists).\n");
                    printf("   probability_file: lines must be formatted as {letter},{nucleotides,3},{weighting},{index}. The weightings do not need to sum to 1. Codon indices must range from 0 to 63.\n");
                    printf("   num_to_subsample: number of 'top' encodings to take, the best encodings are measured by absolute difference from gc_target_ratio.\n" );
                    printf("   number of threads to use for operations, default is 1\n");
                    printf("   gc_target_ratio: ratio GC to AT to target for encodings.\n" );
                    printf("   trials: number of nucleotide sequences to generate for each input sequence. The default is 10,000.\n");
                    printf("   max_line_length: max number of characters allowed in input line. Max line length must not exceed 65,534. Default is 127.\n" );
                    exit(EXIT_SUCCESS);
                default: break;
                }
        }

    // check options
    assert(input_file && probability_file, "Parameter error: need two input files. Type 'codon_sampling -h' for help.\n");
    assert(seq_output_file && ratio_output_file, "Parameter error: output file missing. Type 'codon_sampling -h' for help.\n");
    assert(trials >= 0, "Parameter error: number of trials must be non-negative. Type 'codon_sampling -h' for help.\n");
    // 65534 + 1 to account for the newline char
    assert(max_line_length > 0 && max_line_length < 65536, "Parameter error: max line length must not be 0 or exceed 65,534. Type 'codon_sampling -h' for help.\n" );

    // start
    double begin  = 0.0;

    // seed RNG using hardware source
    xoroshiro::seedrandom();

    // number of digits required in suffix
    size_t digits = (size_t)floor(log10(trials)) + 1;

    // open files
    FILE* fin = fopen(input_file, "r");
    FILE *foutr = fopen( ratio_output_file, "w" );
    if( fin == nullptr )
    {
        fprintf( stderr, "Error opening input file. Verify the file path exists.\n" );
        exit( 1 );
    }
    std::ofstream out_file_seqs;
    auto t = table(probability_file);

    // mapping of amino acid string to index
    char acid_map['Z'] = {0};
    acid_map['A'] = 0;  acid_map['C'] = 1;  acid_map['D'] = 2;  acid_map['E'] = 3;
    acid_map['F'] = 4;  acid_map['G'] = 5;  acid_map['H'] = 6;  acid_map['I'] = 7;
    acid_map['K'] = 8;  acid_map['L'] = 9;  acid_map['M'] = 10; acid_map['N'] = 11;
    acid_map['P'] = 12; acid_map['Q'] = 13; acid_map['R'] = 14; acid_map['S'] = 15;
    acid_map['T'] = 16; acid_map['V'] = 17; acid_map['W'] = 18; acid_map['Y'] = 19;

    // process line by line
    uint64_t lines = count_lines_in_file( input_file );
    uint64_t new_lines = lines;
    uint64_t loop_index = 0;
    uint64_t index = 0;

    char line[ max_line_length ];// = {0};
    char line_copy[ max_line_length ];// = {0};
    std::vector<FileInput> file_data_arr;
    std::vector<std::string> results;

    //omp_set_num_threads( num_threads );

    //replace num_threads

    file_data_arr.reserve( lines );
    results.reserve( lines );

    for( loop_index = 0; loop_index < lines; ++loop_index )
    {
        char* get_line = fgets( line, max_line_length, fin );
        if( get_line == nullptr )
        {
            fprintf( stderr, "Error reading input file.\n" );
            exit( 1 );
        }
        // replace trailing newline with terminator so we can use strlen
        if( strchr( line, '\n' ) == nullptr )
        {
            fprintf( stderr, "Parameter Error: max line length provided has overflowed or is less than provided sequence: %s.\n", strtok( line, "," ) );
            exit( 1 );
        }
        line[strcspn(line, "\n")] = '\0';

        strcpy( line_copy, line );

        // line consists of name,input
        char* input = nullptr;
        char* name = strtok_r( line_copy, ",", &input );


        if( !( strchr( input, 'B' )
               || strchr( input, 'Z' )
               || strchr( input, 'U' )
               || strchr( input, 'O' )
               || strchr( input, 'J' )
            )
          )
            {
                if( valid_line( line, max_line_length ) )
                    {
                        FileInput new_input( name, input );
                        file_data_arr.push_back( new_input );
                    }
                else
                    {
                        report_bad_line( line, loop_index + 1 );
                        --new_lines;
                    }
            }
        else
            {
                printf( "Notice: Skipping oligo with ambiguous code, B,Z U, O or J,  %s.\n",
                        input
                      );
                --new_lines;
            }
        *line = 0;
    }

    lines = new_lines;

    out_file_seqs.open( seq_output_file );

    Encoding **encodings = (Encoding**) malloc( sizeof( Encoding *) * trials );

    for( loop_index = 0; loop_index < lines; ++loop_index )
        {

            FileInput file_data = file_data_arr.at( loop_index );

            file_data.aa_total = file_data.data.length();
            file_data.total_nucleotides = file_data.aa_total * CODON_SIZE;

            for( index = 0; index < file_data.data.length(); ++index )
                {
                    ++file_data.aa_counts[ (uint8_t) acid_map[ (uint8_t) file_data.data[ index ] ] ];
                }


            uint64_t current_aa    = 0;
            uint64_t current_trial = 0;
            uint64_t len = 0;
            uint64_t i = 0;

            len = file_data.data.length();

            Encoding *current = NULL;

            std::mutex mtx;

            //std::mutex mtx;

            unsigned short int threads = num_threads;

            // trials
            //#pragma omp parallel for private( current_trial, current, current_aa ) shared( trials, encodings, len, t ) schedule( static )
            for ( current_trial = 0; current_trial < trials; ++current_trial)
                {

                    std::thread curr_thread = std::thread( [&]{
                        mtx.lock();
                        if (threads == 0) {
                            return;
                        }
                        --threads;
                        ++current_trial;
                        mtx.unlock();
                        current = new Encoding();

                            // keep track of nucleotide and codon ratios
                        current->original = file_data;


                        // calculate result string
                        for ( current_aa = 0; current_aa < len; ++current_aa )
                            {
                                double r = xoroshiro::uniform();
                                const uint16_t aa = file_data.data[current_aa];
                                
                                codon** cod = t[aa];
                                double accum = (*cod)->w;

                                while ( accum < r )
                                    {
                                        accum += (*++cod)->w;
                                    }

                                for ( i = 0; i < 4; ++i )
                                    {
                                        mtx.lock();
                                        current->nucleotides[ i ] += (*cod)->nucleotides[i];
                                        mtx.unlock();
                                    }
                                mtx.lock();
                                ++current->codons[(*cod)->index];
                                mtx.unlock();
                                
                                current->encoding.append( (*cod)->c, CODON_SIZE );

                            }

                        current->calc_gc_ratio();
                        current->total_codons = len;

                        current->gc_dist_abs = fabs( current->gc_ratio - gc_target_ratio );

                        mtx.lock();
                        encodings[ current_trial ] = current;
                        mtx.unlock();

                        threads++;
                    } );

                    curr_thread.join();
                }
            //block until entire for loop is complete
            while ( current_trial != (trials - 1)) {
                
            }

            uint64_t num_encodings = trials;

            // sort the encodings according to gc_ratio
            qsort( encodings, num_encodings, sizeof( Encoding *), encoding_compar );


            std::vector<Encoding *> best_encodings;
            best_encodings.reserve( num_to_subsample );

            for( index = 0; index < num_to_subsample; index++ )
                {
                    best_encodings.push_back( encodings[ index ] );
                }

            std::string out_string;
            for( index = 0; index < num_to_subsample; index++ )
                {
                            std::stringstream digit_str;
                            std::string gc_content;
                            std::string gc_dev;
                            char gc_content_c[ max_line_length ];//
                            char gc_dev_c[ max_line_length ];

                            sprintf( gc_content_c, "%Lf", best_encodings[ index ]->gc_ratio );
                            sprintf( gc_dev_c, "%Lf", best_encodings[ index ]->gc_dist_abs );

                            gc_dev = gc_dev_c;
                            gc_content = gc_content_c;

                            digit_str << std::internal << std::setw( digits ) << std::setfill( '0' ) <<  index + 1;
                            out_string.append( best_encodings[ index ]->original.name );
                            out_string.append( "_" );

                            out_string.append( digit_str.str() );
                            out_string.append( "," );
                            out_string.append( best_encodings[ index ]->original.data );
                            out_string.append( "," );
                            out_string.append( best_encodings[ index ]->encoding );
                            out_string.append( "," );
                            out_string.append( gc_content );
                            out_string.append( "," );
                            out_string.append( gc_dev );
                            out_string.append( "\n" );
                } 

            out_file_seqs << out_string;

            Encoding *current_encoding = NULL;
            char **str_arr = (char**) malloc( sizeof( char *) * num_to_subsample );
            char *new_str = NULL;

            //#pragma omp parallel for shared( str_arr, num_to_subsample) private( new_str, index, current_encoding ) schedule( static )
            for( index = 0; index < num_to_subsample; index++ )
                {

                    current_encoding = best_encodings[ index ];
                    new_str = (char*) malloc( sizeof( char ) * 12 * 88);
                    str_arr[ index ] = new_str;

                    sprintf(new_str, "%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g\n",
                            (double)current_encoding->nucleotides[0]/current_encoding->original.total_nucleotides,
                            (double)current_encoding->nucleotides[1]/current_encoding->original.total_nucleotides,
                            (double)current_encoding->nucleotides[2]/current_encoding->original.total_nucleotides,
                            (double)current_encoding->nucleotides[3]/current_encoding->original.total_nucleotides,
                            (double)current_encoding->original.aa_counts[0]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[1]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[2]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[3]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[4]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[5]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[6]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[7]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[8]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[9]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[10]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[11]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[12]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[13]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[14]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[15]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[16]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[17]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[18]/current_encoding->original.aa_total,
                            (double)current_encoding->original.aa_counts[19]/current_encoding->original.aa_total,
                            (double)current_encoding->codons[0]/current_encoding->total_codons,
                            (double)current_encoding->codons[1]/current_encoding->total_codons,
                            (double)current_encoding->codons[2]/current_encoding->total_codons,
                            (double)current_encoding->codons[3]/current_encoding->total_codons,
                            (double)current_encoding->codons[4]/current_encoding->total_codons,
                            (double)current_encoding->codons[5]/current_encoding->total_codons,
                            (double)current_encoding->codons[6]/current_encoding->total_codons,
                            (double)current_encoding->codons[7]/current_encoding->total_codons,
                            (double)current_encoding->codons[8]/current_encoding->total_codons,
                            (double)current_encoding->codons[9]/current_encoding->total_codons,
                            (double)current_encoding->codons[10]/current_encoding->total_codons,
                            (double)current_encoding->codons[11]/current_encoding->total_codons,
                            (double)current_encoding->codons[12]/current_encoding->total_codons,
                            (double)current_encoding->codons[13]/current_encoding->total_codons,
                            (double)current_encoding->codons[14]/current_encoding->total_codons,
                            (double)current_encoding->codons[15]/current_encoding->total_codons,
                            (double)current_encoding->codons[16]/current_encoding->total_codons,
                            (double)current_encoding->codons[17]/current_encoding->total_codons,
                            (double)current_encoding->codons[18]/current_encoding->total_codons,
                            (double)current_encoding->codons[19]/current_encoding->total_codons,
                            (double)current_encoding->codons[20]/current_encoding->total_codons,
                            (double)current_encoding->codons[21]/current_encoding->total_codons,
                            (double)current_encoding->codons[22]/current_encoding->total_codons,
                            (double)current_encoding->codons[23]/current_encoding->total_codons,
                            (double)current_encoding->codons[24]/current_encoding->total_codons,
                            (double)current_encoding->codons[25]/current_encoding->total_codons,
                            (double)current_encoding->codons[26]/current_encoding->total_codons,
                            (double)current_encoding->codons[27]/current_encoding->total_codons,
                            (double)current_encoding->codons[28]/current_encoding->total_codons,
                            (double)current_encoding->codons[29]/current_encoding->total_codons,
                            (double)current_encoding->codons[30]/current_encoding->total_codons,
                            (double)current_encoding->codons[31]/current_encoding->total_codons,
                            (double)current_encoding->codons[32]/current_encoding->total_codons,
                            (double)current_encoding->codons[33]/current_encoding->total_codons,
                            (double)current_encoding->codons[34]/current_encoding->total_codons,
                            (double)current_encoding->codons[35]/current_encoding->total_codons,
                            (double)current_encoding->codons[36]/current_encoding->total_codons,
                            (double)current_encoding->codons[37]/current_encoding->total_codons,
                            (double)current_encoding->codons[38]/current_encoding->total_codons,
                            (double)current_encoding->codons[39]/current_encoding->total_codons,
                            (double)current_encoding->codons[40]/current_encoding->total_codons,
                            (double)current_encoding->codons[41]/current_encoding->total_codons,
                            (double)current_encoding->codons[42]/current_encoding->total_codons,
                            (double)current_encoding->codons[43]/current_encoding->total_codons,
                            (double)current_encoding->codons[44]/current_encoding->total_codons,
                            (double)current_encoding->codons[45]/current_encoding->total_codons,
                            (double)current_encoding->codons[46]/current_encoding->total_codons,
                            (double)current_encoding->codons[47]/current_encoding->total_codons,
                            (double)current_encoding->codons[48]/current_encoding->total_codons,
                            (double)current_encoding->codons[49]/current_encoding->total_codons,
                            (double)current_encoding->codons[50]/current_encoding->total_codons,
                            (double)current_encoding->codons[51]/current_encoding->total_codons,
                            (double)current_encoding->codons[52]/current_encoding->total_codons,
                            (double)current_encoding->codons[53]/current_encoding->total_codons,
                            (double)current_encoding->codons[54]/current_encoding->total_codons,
                            (double)current_encoding->codons[55]/current_encoding->total_codons,
                            (double)current_encoding->codons[56]/current_encoding->total_codons,
                            (double)current_encoding->codons[57]/current_encoding->total_codons,
                            (double)current_encoding->codons[58]/current_encoding->total_codons,
                            (double)current_encoding->codons[59]/current_encoding->total_codons,
                            (double)current_encoding->codons[60]/current_encoding->total_codons,
                            (double)current_encoding->codons[61]/current_encoding->total_codons,
                            (double)current_encoding->codons[62]/current_encoding->total_codons,
                            (double)current_encoding->codons[63]/current_encoding->total_codons);
                }

            for( index = 0; index < num_to_subsample; ++index )
                {
                            fprintf( foutr, "%s", str_arr[ index ] );
                            free( str_arr[ index ] );
                }
            free( str_arr );

            for( index = 0; index < trials; index++ )
                {
                    delete encodings[ index ];
                }
        }

    free( encodings );

    out_file_seqs.close();

    fclose(fin);

    printf("Processed %u lines x %lu trials, time elapsed %f s\n", lines, trials, begin - begin );

    return EXIT_SUCCESS;
}

uint64_t count_lines_in_file( const char *filename )
{
    FILE *open_file = NULL;
    int c = 0;
    uint64_t count = 0;

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

int valid_input( char *input )
{
    int cur_char = 0;

    if( input )
        {
            while( input[ cur_char ] )
                {
                    if( !( input[ cur_char ] >= 'A' &&
                           input[ cur_char ] <= 'Z'
                         )
                      )
                        {
                            // input contains invalid char
                            return 0;
                        }
                    ++cur_char;
                }
            // input is valid
            return 1;
        }
    // input is null
    return 0;
}

int valid_line( char *input_line, uint16_t max_line_length )
{
    char copy_str[ max_line_length ];
    strcpy( copy_str, input_line );

    // line consists of name,input
    char* input = nullptr;
    char* name = strtok_r( copy_str, ",", &input );

    if( name && valid_input( input ) )
        {
            return 1;
        }
    return 0;
}

void report_bad_line( char *line, int lineno )
{
    printf( "Warning: Line %d: %s is invalid and will be skipped.\n", lineno, line );
}
