#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "protein_oligo_library.h"
#include "dynamic_string.h"
#include "array_list.h"

#define LINE_SIZE 512
#define DASH_CHAR '-'
#define SPACE ' '
#define DATA_NOT_FOUND -99

// ============== Local Function Prototypes ==================== // 

char get_first_char_in_functional_group( char in_char );
char get_corresponding_char( char in_char );
void xmer_first_functional_group( char* in_string, int str_len  );
int count_letters( char* in_str );
void get_blosum_distances( hash_table_t* blosum_distance, FILE* blosum_file, int num_rows );
void get_alpha_chars( char* dest, char* source, int num_chars );
void get_ints_from_string( int* dest, char* src, int num_rows );
int get_blosum_dist( blosum_data_t* in_data, char first, char second );
void write_xmer( char* to_write, char* in_seq, int xmer_index, int window_size, int step_size );



// ============================================================= // 

int num_digits_in_int( int input )
{
    char int_as_str[ LINE_SIZE ];
    return sprintf( int_as_str, "%d", input );
}

int count_letters( char* in_str )
{
    int char_count = 0;
    int index = 0;

    while( in_str[ index ] )
        {
            if( in_str[ index ] >= 'A' && in_str[ index ] <= 'Z' )
                {
                    char_count++;
                }
            index++;
        }

    return char_count;
}

int get_blosum_dist( blosum_data_t* in_data, char first, char second )
{
    int index = 0;
    int return_dist = DATA_NOT_FOUND;
    int len = strlen( in_data->letter_data );
    int* first_distances = ht_find( in_data->blosum_table, &first );

    if( first_distances != NULL )
        {
            while( in_data->letter_data[ index ] != second &&
                   index < len
                 )
                {
                    index++;
                }

            return_dist = first_distances[ index ];
        }

    return return_dist;
}

void get_alpha_chars( char* dest, char* source, int num_chars )
{
    int char_count = 0;
    int index = 0;

    while( index < num_chars )
        {
            if( source[ char_count ] >= 'A' && source[ char_count ] <= 'Z' )
                {
                    dest[ index ] = source[ char_count ];
                    index++;
                }
            char_count++;
        }
    dest[ index ] = '\0';
}

void get_ints_from_string( int* dest, char* src, int num_rows )
{
    int index = 0;
    char* token = NULL;

    // eat leading letter
    token = strtok( src, " " );
    token = strtok( NULL, " " );

    while( token && index < num_rows )
        {
            dest[ index ] = atoi( token );
            token = strtok( NULL, " " );
            index++;
        }
}

void get_blosum_distances( hash_table_t* blosum_distance, FILE* blosum_file, int num_rows )
{
    char* found_line = NULL;
    char current_line[ LINE_SIZE ] ;
    int* current_distance_data = NULL;
    char key[ 2 ];

    found_line = fgets( current_line, LINE_SIZE, blosum_file );

    while( found_line )
        {
            current_distance_data = malloc( sizeof( int ) * num_rows );

            strncpy( key, current_line, 1 );
            key[ 1 ] = '\0';

            get_ints_from_string( current_distance_data, current_line, num_rows );
            ht_add( blosum_distance, key, current_distance_data );

            found_line = fgets( current_line, LINE_SIZE, blosum_file );
        }
}

void read_sequences( FILE* file_to_read, sequence_t** in_sequence )
{
    int has_line;
    int index = 0;

    dynamic_string_t* line = (dynamic_string_t*) malloc( sizeof( dynamic_string_t ) );
    dynamic_string_t* sequence = malloc( sizeof( dynamic_string_t ) );

    has_line = get_a_line( file_to_read, line );
    while( has_line )
        {

            if( line->data[ 0 ] == '>' )
                {
                    *( in_sequence + index ) = malloc( sizeof( sequence_t ) );
                    sequence = malloc( sizeof( dynamic_string_t ) );
                    ds_init( sequence );

                    in_sequence[ index ]->name = line->data;
                    in_sequence[ index ]->sequence = sequence; 
                    index++;
                }
            else
                {
                    ds_add( sequence, line->data );
                }
            ds_init( line );
            has_line = get_a_line( file_to_read, line );
        }

    ds_clear( line );
}

void write_fastas( sequence_t** in_seqs, int num_seqs, char* output )
{
    FILE* out_file = fopen( output, "w+" );
    int index;

    if( !out_file )
        {
            printf( "Unable to open file %s for output." , output );
        }

    for( index = 0; index < num_seqs; index++ )
        {
            fprintf( out_file, "%s\n", in_seqs[ index ]->name );

            fprintf( out_file, "%s\n", in_seqs[ index ]->sequence->data );
        }

    fclose( out_file );
}


int count_seqs_in_file( FILE* data_file )
{

    int counter = 0;
    int current_char = 0;

    if( !data_file )
        {
            return -1; 
        }

    rewind( data_file );

    current_char = fgetc( data_file );
    while( current_char != EOF )
        {
            if( (char) current_char == '>' )
                {
                    counter++;
                }
            current_char =  fgetc( data_file );
        }
    rewind( data_file );
    return counter;
}

void write_xmer( char* to_write, char* in_seq, int xmer_index, int window_size, int step_size )
{
    int index = 0;

    for( index = 0; index < window_size; index++ )
        {
            to_write[ index ] =
                in_seq[ ( xmer_index * step_size ) + index ];
        }

    to_write[ index ] = '\0';
}
int get_a_line( FILE* stream, dynamic_string_t* to_read )
{
    char current_char[ 256 ] ;

    ds_init( to_read );
    if( fgets( current_char, 256, stream ) ) 
        {
            ds_add( to_read, current_char );
            return true;
        }
    return false;

}


int count_char_in_string( char* string_in, char to_find )
{
    int count = 0;
    int index = 0;

    while( *(string_in + index ) )
        {
            if( *( string_in + index ) == to_find )
                {
                    count++;
                }
            index++;
        }

    return count;
}

int char_in_string( char* string_in, char to_find )
{
    int index = 0;
    while( *( string_in + index ) )
        {
            if ( *(string_in + index ) == to_find )
                {
                    return true;
                }
            index++;
        }
    return false;
}

float percent_char_in_string( char* string_in, char test_char )
{
    return ( (float) count_char_in_string( string_in, test_char ) /
             string_length( string_in ) ) * 100 ;
}

void remove_char_from_string( char* test_string, char to_remove )
{
    int index;
    int valid_index = 0;

    for( index = 0; test_string[ index ]; index++ )
        {
            if( test_string[ index ] != to_remove )
                {
                    test_string[ valid_index ] = test_string[ index ];
                    valid_index++;
                }
        }
    *( test_string + valid_index ) = '\0';
}

int is_valid_sequence( char* sequence, int min_length, float percent_valid )
{
    if( !char_in_string( sequence, 'X' ) )
        {
            if( !min_length )
                {
                    return percent_char_in_string( sequence, DASH_CHAR )
                           < ( 100 - percent_valid );
                }
            return count_char_in_string( sequence, DASH_CHAR )
                <= ( string_length( sequence ) - min_length );
        }
   return 0;
}


int calc_num_subseqs( int length, int window_size )
{
    int return_val = 0;
    if( length >= window_size )
        {
            return_val = length - window_size + 1;
        }

    return return_val;
}

void append_suffix( char* result, char* in_name, int start, int end )
{
    sprintf( result, "%s_%d_%d", in_name, start, end );
}

hash_table_t* subset_lists( hash_table_t* in_hash,
                            char* in_seq,
                            int window_size, int step_size,
                            blosum_data_t* blosum_data,
                            int blosum_cutoff,
                            int permute
                          )
{
    int outer_index;
    int num_subsets = calc_num_subseqs( strlen( in_seq ), window_size );

    uint32_t permute_index;

    char *current_xmer;

    array_list_t *current_xmer_permutations;

    for( outer_index = 0; outer_index < num_subsets; outer_index++ )
        {

            current_xmer = malloc( window_size + 1 );


            write_xmer( current_xmer, in_seq, outer_index, window_size, step_size );

            if( permute || blosum_data )
                {
                    current_xmer_permutations = malloc( sizeof( array_list_t ) );
                    ar_init( current_xmer_permutations );
                    permute_xmer_functional_groups( current_xmer, current_xmer_permutations,
                                                    blosum_data, blosum_cutoff
                                                  );

                    for(  permute_index = 0; permute_index < current_xmer_permutations->size; permute_index++ )
                        {
                            current_xmer = ar_get( current_xmer_permutations, permute_index );
                            ht_add( in_hash, current_xmer, NULL );

                            free( current_xmer );
                        }

                    ar_clear( current_xmer_permutations );

                }

            // update the entry at this location
            ht_add( in_hash, current_xmer, NULL );

        }
    return in_hash;
}

hash_table_t* create_xmers_with_locs( hash_table_t* in_hash, char* in_name,
                                      char* in_seq,
                                      int window_size, int step_size )
{
    int outer_index;
    int num_subsets = calc_num_subseqs( strlen( in_seq ), window_size );

    subset_data_t xmer_bounds[ num_subsets ];
    subset_data_t current_xmer_data;

    array_list_t* xmer_locations;

    char *current_xmer;
    char* name_with_bounds;

    for( outer_index = 0; outer_index < num_subsets; outer_index++ )
        {

            current_xmer_data = xmer_bounds[ outer_index ];

            current_xmer = malloc( window_size + 1 );
            current_xmer_data.start = ( outer_index * step_size );

            write_xmer( current_xmer, in_seq, outer_index, window_size, step_size );

            current_xmer_data.end = ( outer_index * step_size ) + window_size;

            name_with_bounds = malloc( strlen( in_name ) +
                                       num_digits_in_int( current_xmer_data.start ) +
                                       num_digits_in_int( current_xmer_data.end ) +
                                       2+ 1
                                      );


            append_suffix( name_with_bounds, in_name, current_xmer_data.start, current_xmer_data.end );

            if( in_hash != NULL &&
			    !char_in_string( current_xmer, 'X' )
              ) 
            {
                xmer_locations = ( array_list_t* ) ht_find( in_hash, current_xmer );
                if( xmer_locations != NULL )
                    {
                        // update the entry at this location
                        ar_add( xmer_locations, name_with_bounds );
                    }
                else
                    {
                        // create the entry at this location
                        xmer_locations = malloc( sizeof( array_list_t ) );

                        ar_init( xmer_locations );
                        ar_add( xmer_locations, name_with_bounds );
 
                        ht_add( in_hash, current_xmer, xmer_locations );
                    }
            }
            else
                {
                    free( current_xmer );
                }
        }
    return in_hash;
}




set_t* component_xmer_locs( char* in_ymer_name, char* in_ymer,
                            set_t* out_ymer,
                            hash_table_t* in_xmer_table,
                            int window_size, int step_size,
                            blosum_data_t* blosum_data,
                            int blosum_cutoff, 
                            int permute
                          )
{
    int num_xmers = ( window_size - step_size ) + 1;
    uint32_t inner_index;
    uint32_t index;
    hash_table_t* subset_xmers = NULL;
    HT_Entry* subset_xmer_items = NULL;
    array_list_t* found_data = NULL;

    uint32_t size;

    subset_xmers = malloc( sizeof( hash_table_t ) );

    ht_init( subset_xmers, num_xmers );

    create_xmers_with_locs( subset_xmers, in_ymer_name, in_ymer,
                            window_size, step_size );

    subset_xmer_items = ht_get_items( subset_xmers );

    size = subset_xmers->size;

    if( permute || blosum_data )
        {
            for( index = 0; index < size; index++ )
                {

                    found_data = malloc( sizeof( array_list_t ) );
                    ar_init( found_data );

                    permute_xmer_functional_groups( subset_xmer_items[ index ].key, found_data,
                                                    blosum_data, blosum_cutoff
                                                  );
                    for( inner_index = 0; inner_index < found_data->size; inner_index++ )
                        {
                            ht_add( subset_xmers, ar_get( found_data, inner_index ), NULL );
                        }


                    for( inner_index = 0; inner_index < found_data->size; inner_index++ )
                        {
                            free( ar_get( found_data, inner_index ) );
                        }
                    ar_clear( found_data );
                }

        }

   free( subset_xmer_items );
   subset_xmer_items = ht_get_items( subset_xmers );

    for( index = 0; index < subset_xmers->size; index++ )
        {
            found_data = (array_list_t*) ht_find( in_xmer_table, subset_xmer_items[ index ].key );
            if( found_data != NULL )
                {
                    set_add_all( out_ymer, (char**) found_data->array_data, found_data->size );

                    if( subset_xmer_items[ index ].value )
                        {
                            ar_clear( subset_xmer_items[ index ].value );
                        }
                }
        }

    free( subset_xmer_items );
    ht_clear( subset_xmers );
    free( subset_xmers );
    return out_ymer;
}


void *malloc_track( array_list_t* data, int num_bytes )
{
    void* tracked_ptr = malloc( num_bytes );
    ar_add( data, tracked_ptr );

    return tracked_ptr;
}

void free_data( array_list_t* in_data )
{
    unsigned int index;
    for( index = 0; index < in_data->size; index++ )
        {
            free( in_data->array_data[ index ] );
        }
}

void permute_xmer_functional_groups( char* str_to_change,
                                     array_list_t* permutations,
                                     blosum_data_t* blosum_data,
                                     int blosum_cutoff
                                   )
{
    int length = strlen( str_to_change );
    int index;
    size_t inner_index;
    int blosum_dist;

    char *copied_string = NULL;
    char original_char = '\0';
    char copy_string[ length + 1 ];

    char different_char = '\0';

    copied_string = malloc( sizeof( char ) * length + 1 );
    for( index = 0; index < length; index++ )
        {
            memset( copy_string, '\0', sizeof( copy_string ) );
            strcpy( copy_string, str_to_change );
            original_char = copy_string[ index ];

            different_char = get_first_char_in_functional_group( copy_string[ index ] );

            if( blosum_data )
                {
                    for( inner_index = 0; inner_index < strlen( blosum_data->letter_data ); inner_index++ )
                        {
                            different_char = blosum_data->letter_data[ inner_index ];
                            blosum_dist = get_blosum_dist( blosum_data, original_char, different_char );
                            if( blosum_dist >= blosum_cutoff )
                                {
                                    copy_string[ index ] = different_char;
                                    copied_string = malloc( sizeof( char ) * length + 1 );
                                    strcpy( copied_string, copy_string );
                                    ar_add( permutations, copied_string );
                                }
                        }
                }
            else
                {
                    while( different_char )
                        {
                            copy_string[ index ] = different_char;
                            copied_string = malloc( sizeof( char ) * length + 1 );
                            strcpy( copied_string, copy_string );
                            ar_add( permutations, copied_string );

                            different_char = get_corresponding_char( different_char );
                        }

                }
        }
}

void xmer_first_functional_group( char* in_string, int str_len  )
{
    int index;
    
    for( index = 0; index < str_len; index++ )
        {
            in_string[ index ] = get_first_char_in_functional_group( in_string[ index ] );
        }
}

char get_corresponding_char( char in_char )
{
    switch( in_char )
        {
        case 'E':
            return 'D';

        case 'H':
            return 'K';
        case 'K':
            return 'R';

        case 'C':
            return 'T';
        case 'T':
            return 'S';
        case 'S':
            return 'N';
        case 'N':
            return 'Q';

        case 'F':
            return 'Y';
        case 'Y':
            return 'W';

        case 'A':
            return 'V';
        case 'V':
            return 'M';
        case 'M':
            return 'L';
        case 'L':
            return 'I';
            
        default:
            return '\0';
        }
}

char get_first_char_in_functional_group( char in_char )
{
    switch( in_char )
        {
        case 'R': // FALLTHROUGH INTENTIONAL
        case 'K':
            return 'H';
        case 'E':
            return 'D';
        case 'T': // FALLTHROUGH INTENTIONAL
        case 'S':
        case 'Q':
        case 'N':
            return 'C';
        case 'Y': // FALLTHROUGH INTENTIONAL
        case 'W':
            return 'F';
        case 'V': // FALLTHROUGH INTENTIONAL
        case 'M':
        case 'L':
        case 'I':
            return 'A';
        default:
            return '\0';
        }
}


blosum_data_t* parse_blosum_file( FILE* file_name )
{
    FILE* blosum = file_name;

    int row_count = 0;

    char* letter_data = NULL;
    char* current_line = malloc( LINE_SIZE );

    hash_table_t* blosum_distances = NULL;
    blosum_data_t* blosum_data = NULL;

    if( blosum )
        {
            current_line = fgets( current_line, LINE_SIZE, blosum );

            // consume any preceeding characters
            while( current_line && current_line[ 0 ] == '#' )
                {
                    current_line = fgets( current_line, LINE_SIZE, blosum );
                }
            if( current_line )
                {
                    row_count = count_letters( current_line );

                    letter_data = malloc( sizeof( char ) * row_count + 1 );
                    blosum_data = malloc( sizeof( blosum_data_t ) );

                    blosum_distances = malloc( sizeof( hash_table_t ) );
                    ht_init( blosum_distances, row_count );

                    get_alpha_chars( letter_data, current_line, row_count );
                    get_blosum_distances( blosum_distances, blosum, row_count );

                    blosum_data->letter_data = letter_data;
                    blosum_data->blosum_table = blosum_distances;
                    
                }

                
        }
    free( current_line );

    return blosum_data;
}



