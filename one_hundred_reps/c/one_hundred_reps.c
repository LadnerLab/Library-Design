#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>

#include "protein_oligo_library.h"
#include "array_list.h"
#include "hash_table.h"

const int NUM_THREADS   = 4;
const int MAX_NUM_CHARS = 256;
const char *ARGS        = "f:o:h::n:m:";
const int SMALL_TABLE_SIZE = 1000;

int seq_compare( const void *a, const void *b );
void write_map( hash_table_t *table, char *filename );

int main( int argc, char **argv )
{
    int num_threads = NUM_THREADS;
    int option      = 0;
    int index       = 0;
    int inner_index = 0;
    int outer_index = 0;
    int num_seqs    = 0;

    char in_file[ MAX_NUM_CHARS ];
    char outfile_name[ MAX_NUM_CHARS ];

    FILE *file = NULL;

    bool found = false;

    sequence_t **in_seqs;

    array_list_t *out_seqs = NULL;
    hash_table_t *new_list = NULL;

    sequence_t **intermed_seqs = NULL;
    hash_table_t *map_table    = NULL;
    char map_file[ MAX_NUM_CHARS ];
    int TABLE_SIZE = 0;

    hash_table_t *other_table = NULL;
    HT_Entry **other_entries = NULL;
    uint64_t innermost_index = 0;
    bool map_file_included = false;

    out_seqs = malloc( sizeof( array_list_t ) );

    ar_init( out_seqs );

    while( ( option = getopt( argc, argv, ARGS ) ) != -1 )
        {
            switch( option )
                {
                    case 'f': // infile_name
                        strcpy( in_file, optarg );
                        break;
                    case 'o': // outfile_name
                        strcpy( outfile_name, optarg );
                        break;
                    case 'n': // num_threads
                        num_threads = atoi( optarg );
                        break;
                    case 'm': // map_file
                        strcpy( map_file, optarg );
                        map_file_included = true;
                        break;
                    case 'h':
                        printf( "USAGE:\nexecutable: 'one_hundred_reps'\n"
                                "-f                            Name of fasta file to parse and collapse.\n");
                        printf( "-o                            Name of the output fasta file.\n" );
                        printf( "-n                            Number of threads to use. Default is 4.\n" );
                        printf( "-m                            Optional, name of file to write out map to, this map will\n"
                                "                              contain a mapping of which sequences were collapsed under\n"
                                "                              which sequences.\n" );
                        printf( "-h                            Provides a help message.");
                        return EXIT_SUCCESS;
                    
                    default:
                        printf( "Incorrect argument supplied!\n" );

                        return EXIT_FAILURE;
                    
                }
        }

    file = fopen( in_file, "r" );



    omp_set_num_threads( num_threads );

    if( !file )
        {
            printf( "ERROR: %s could not be found, program will exit...\n", in_file );
            return EXIT_FAILURE;
        }

    num_seqs = count_seqs_in_file( file );
    TABLE_SIZE = num_seqs;

    printf( "Num Seqs: %d\n", num_seqs );


    in_seqs  = malloc( sizeof( sequence_t * ) * num_seqs );
    intermed_seqs = calloc( num_seqs, sizeof( sequence_t *) );

    read_sequences( file, in_seqs );

    fclose( file );

    sequence_t *copy_seqs = malloc( sizeof( sequence_t ) * num_seqs );

    for( index = 0; index < num_seqs; index++ )
        {
            copy_seqs[ index ] = *in_seqs[ index ];
        }

    qsort( copy_seqs, num_seqs, sizeof( sequence_t ), seq_compare );

    for( index = 0; index < num_seqs; index++ )
        {
            in_seqs[ index ] = &copy_seqs[ index ];
        }


    if( map_file_included )
        {
            // map table populated by in_seqs/input fasta
            map_table = malloc( sizeof( hash_table_t ) );
            ht_init( map_table, TABLE_SIZE );

            for( index = 0; index < num_seqs; index++ )
                {
                    new_list = malloc( sizeof( hash_table_t ) );
                    ht_init( new_list, SMALL_TABLE_SIZE );
                    ht_add( new_list, in_seqs[ index ]->name, in_seqs[ index ] );
                    ht_add( map_table, in_seqs[ index ]->name, new_list );
                }
        }
    // To Do: Rewrite this to use 3 maps: original in_seq map (seq name : sequence), new no dup in_seq map (seq name : sequence), coverage map (cover seq name : list of covered seq names)
    #pragma omp parallel for private( other_table, other_entries, innermost_index, outer_index, inner_index, found, new_list ) \
                             shared( intermed_seqs, in_seqs, num_seqs, map_table ) schedule( dynamic )
    for( outer_index = 0; outer_index < num_seqs; outer_index++ )
        {
            found = false;
            for( inner_index = num_seqs - 1; inner_index > outer_index; inner_index-- )
                {
                    // if pointer to occurrence of second seq in first seq returned
                    if( strstr( in_seqs[ inner_index]->sequence->data,
                                in_seqs[ outer_index ]->sequence->data
                              )
                      )
                        {
                            found = true;

                            if( map_table )
                                {
                                    #pragma omp critical
                                    {
                                        // if not collapsed
                                        if( in_seqs[ outer_index ]->collapsed == 0)
                                            {
                                                // set collapsed
                                                in_seqs[ outer_index ]->collapsed = 1;
                                                
                                                other_table = ht_delete( map_table, in_seqs[ outer_index ]->name );

                                                new_list = ht_find( map_table, in_seqs[ inner_index ]->name );
                                                other_entries = ht_get_items( other_table );
                                                for( innermost_index = 0; innermost_index < other_table->size; innermost_index++ )
                                                    {
                                                        ht_add( new_list, other_entries[ innermost_index ]->key, other_entries[ innermost_index ]->value );
                                                    }

                                                free( other_entries );
                                            }
                                    }
                                }
                            break;
                        }

                }
            if( !found )
                {
                    intermed_seqs[ outer_index ] = in_seqs[ outer_index ];
                }
        }

    for( index = 0; index < num_seqs; index++ )
        {
            if( intermed_seqs[ index ] && intermed_seqs[ index ]->collapsed == 0 )
                {
                    ar_add( out_seqs, intermed_seqs[ index ] );
                }
        }
    printf( "Output seqs: %d\n", out_seqs->size );

    if( map_file_included )
        {
            write_map( map_table, map_file );
        }

    write_fastas( (sequence_t**)out_seqs->array_data, out_seqs->size, outfile_name );

    return EXIT_SUCCESS;
}

int seq_compare( const void *a, const void *b )
{
    const sequence_t *first  = a;
    const sequence_t *second = b;

    return first->sequence->size - second->sequence->size;
}

void write_map( hash_table_t *table, char *filename )
{
    HT_Entry **ht_items        = NULL;
    hash_table_t *current_arr = NULL;
    HT_Entry **current_arr_items = NULL;
    FILE *out_file            = NULL;

    uint32_t index       = 0;
    uint32_t inner_index = 0;

    // sequence_t *current_data = NULL;

    out_file = fopen( filename, "w" );

    if( out_file )
        {
            ht_items = ht_get_items( table );
            for( index = 0; index < table->size; index++ )
                {
                    current_arr = ht_items[ index ]->value; // covered seq names

                    current_arr_items = ht_get_items( current_arr );
                    fprintf( out_file, "%s\t", ++( ht_items[ index ]->key ) ); // output fasta seq name to first column

                    if( current_arr->size > 1 )
                        {
                            for( inner_index = 0; inner_index < current_arr->size; inner_index++ ) // 3rd index, 18136 seq name segfaults. This is the sequence that gets left out and a duplicate of the cover seq name takes its place normally
                                {
                                    if( strstr( current_arr_items[ inner_index ]->key, ht_items[ index ]->key ) == NULL )
                                        {
                                            // current_data = current_arr_items[ inner_index ]->value;
                                            fprintf( out_file, "%s\t", ++(current_arr_items[ inner_index ]->key) );
                                        }
                                }

                        }
                    fprintf( out_file, "\n" );

                    free( current_arr_items );
                }
            free( ht_items );
            fclose( out_file );
        }
}
