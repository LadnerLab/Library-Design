#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "set.h"
#include "hash_table.h"

#define DEFAULT_SIZE 1000

void set_init( set_t* to_init )
{
    to_init->data = malloc( sizeof( hash_table_t ) );
    ht_init( to_init->data, DEFAULT_SIZE );
}


void set_add( set_t* set_to_add, char* add_data )
{
   ht_add( set_to_add->data, add_data, 0 );
}

int set_remove( set_t* set_to_remove, char* remove_data )
{
    void* return_result = ht_delete( set_to_remove->data, remove_data );

    free( return_result );

    return return_result == NULL;
}

void set_clear( set_t* set_to_clear )
{
    ht_clear( set_to_clear->data );
    free( set_to_clear->data );
    free( set_to_clear );
}

void set_difference( set_t* first, set_t* second )
{
    uint32_t index;
    uint32_t max = first->data->size;

    HT_Entry **found_data;

    found_data = ht_get_items( first->data );

    for( index = 0; index < ( max ); index++ )
        {
                   
            if( find_item( second->data, found_data[ index ]->key ) )
                {
                    free( ht_delete( first->data, found_data[ index ]->key ) );
                }

        }
    free( found_data );
}

void set_add_all( set_t* dest, char** in_array, int num_elements )
{
    int index;

    for( index = 0; index < num_elements; index++ )
        {
            set_add( dest, in_array[ index ] );
        }
}
