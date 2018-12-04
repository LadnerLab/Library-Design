#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>

#include "protein_oligo_library.h"
#include "array_list.h"

const int NUM_THREADS   = 4;
const int MAX_NUM_CHARS = 256;
const char *ARGS        = "f:n:";

int seq_compare( const void *a, const void *b );

int main( int argc, char **argv )
{
    int num_threads = NUM_THREADS;
    int option      = 0;
    int index       = 0;
    int inner_index = 0;
    int outer_index = 0;
    int num_seqs    = 0;

    char in_file[ MAX_NUM_CHARS ];
    char out_file[ MAX_NUM_CHARS ];

    FILE *file = NULL;

    bool found = false;

    sequence_t **in_seqs;

    array_list_t *out_seqs = NULL;

    out_seqs = malloc( sizeof( array_list_t ) );

    ar_init( out_seqs );

    while( ( option = getopt( argc, argv, ARGS ) ) != -1 )
        {
            switch( option )
                {
                    case 'f':
                        strcpy( in_file, optarg );
                        strcpy( out_file, optarg );
                        strcat( out_file, "_out" );
                        break;
                    case 'n':
                        num_threads = atoi( optarg );
                        break;
                    
                }
        }

    file = fopen( in_file, "r" );

    if( !file )
        {
            printf( "ERROR: %s could not be found, program will exit...\n", in_file );
            return EXIT_FAILURE;
        }

    num_seqs = count_seqs_in_file( file );

    printf( "Num Seqs: %d\n", num_seqs );

    in_seqs  = malloc( sizeof( sequence_t * ) * num_seqs );

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

    for( outer_index = 1; outer_index < num_seqs; outer_index ++ )
        {
            found = false;
            for( inner_index = 0; inner_index < outer_index; inner_index++ )
                {
                    if( strstr( in_seqs[ outer_index ]->sequence->data,
                                in_seqs[ inner_index ]->sequence->data
                              )
                      )
                        {
                            found = true;
                            break;
                        }
                }
            if( !found )
                {
                    ar_add( out_seqs, in_seqs[ outer_index ] );
                }
                 
        }

    printf( "Output seqs: %d\n", out_seqs->size );
    for( index = 0; index < num_seqs; index++ )
        {
            ds_clear( in_seqs[ index ]->sequence );
            free( in_seqs[ index ]->sequence );
        }

    free( in_seqs );
    free( out_seqs );

    return EXIT_SUCCESS;
}

int seq_compare( const void *a, const void *b )
{
    const sequence_t *first  = a;
    const sequence_t *second = b;

    return first->sequence->size - second->sequence->size;
}
