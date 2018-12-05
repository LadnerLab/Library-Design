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
const char *ARGS        = "f:n:m:";

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
    array_list_t *new_list = NULL;

    sequence_t **intermed_seqs = NULL;
    hash_table_t *map_table    = NULL;
    HT_Entry *items = NULL;
    char map_file[ MAX_NUM_CHARS ];
    int TABLE_SIZE = 0;

    bool map_file_included = false;

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
                    case 'm':
                        map_file_included = true;
                        break;
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

    if( map_file_included )
        {
            map_table = malloc( sizeof( hash_table_t ) );
            ht_init( map_table, TABLE_SIZE );
        }

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

    #pragma omp parallel for private( outer_index, inner_index, found ) shared( intermed_seqs, in_seqs, num_seqs ) schedule( dynamic )
    for( outer_index = 0; outer_index < num_seqs; outer_index++ )
        {
            found = false;
            for( inner_index = outer_index + 1; inner_index < num_seqs; inner_index++ )
                {
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
                                        in_seqs[ outer_index ]->collapsed = 1;

                                        if( in_seqs[ inner_index ]->collapsed == 0 )
                                            {
                                                new_list = ht_find( map_table, in_seqs[ inner_index ]->name );
                                                if( !new_list )
                                                    {
                                                        new_list = malloc( sizeof( array_list_t ) );
                                                        ar_init( new_list );
                                                        ht_add( map_table, in_seqs[ inner_index ]->name, new_list );

                                                    }
                                                ar_add( new_list, in_seqs[ outer_index ] );
                                            }

                                    }
                                }
                            break;
                        }

                }
            if( !found )
                intermed_seqs[ outer_index ] = in_seqs[ outer_index ];
        }

    for( index = 0; index < num_seqs; index++ )
        {
            if( intermed_seqs[ index ] )
                {
                    ar_add( out_seqs, intermed_seqs[ index ] );
                }
        }
    printf( "Output seqs: %d\n", out_seqs->size );

    write_fastas( (sequence_t**)out_seqs->array_data, out_seqs->size, out_file );
    for( index = 0; index < num_seqs; index++ )
        {
            ds_clear( in_seqs[ index ]->sequence );
            free( in_seqs[ index ]->sequence );
        }

    free( in_seqs );
    free( out_seqs );

    if( map_file_included )
        {
            items = ht_get_items( map_table );
            uint32_t table_index = 0;

            for( table_index = 0; table_index < map_table->size; table_index++ )
                {
                    ar_clear( items[ table_index ].value );
                }
            ht_clear( map_table );
            free( items );
            free( map_table );
        }

    return EXIT_SUCCESS;
}

int seq_compare( const void *a, const void *b )
{
    const sequence_t *first  = a;
    const sequence_t *second = b;

    return first->sequence->size - second->sequence->size;
}
