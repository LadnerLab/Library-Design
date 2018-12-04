#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <string.h>

#include "protein_oligo_library.h"
#include "array_list.h"

const int NUM_THREADS   = 4;
const int MAX_NUM_CHARS = 256;
const char *ARGS        = "f:n:";

int main( int argc, char **argv )
{
    int num_threads = NUM_THREADS;
    int option      = 0;
    int index       = 0;

    char in_file[ MAX_NUM_CHARS ];
    char out_file[ MAX_NUM_CHARS ];


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

    return EXIT_SUCCESS;
}
