#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "hash_table.h"

#define HASH_NUMBER 3187
#define ADDITIONAL_SPACE 256

// local method for calculating exponents
int int_to_pow( int base, int exponent )
{
    if( exponent == 0 )
        {
            return 1;
        }
    return base * int_to_pow( base, exponent - 1 );
}


void ht_init( hash_table_t* table, int size )
{
    int index;

    table->table_data = malloc( sizeof( HT_Entry* ) * ( size + ADDITIONAL_SPACE ) ); 
    table->capacity = size;
    table->size = 0;

    for( index = 0; index < size; index++ )
        {
            table->table_data[ index ] = NULL;
        }

}

void ht_clear( hash_table_t* table )
{
    uint32_t index;
    HT_Entry* current_node;
    for( index = 0; index < table->capacity; index++ )
        {
            current_node = table->table_data[ index ];
            if( current_node != NULL )
                {
                    while( current_node->next != NULL )
                        {
                            current_node = current_node->next;
                        }
                    while( current_node->prev != NULL )
                        {
                            current_node = current_node->prev;
                            free( current_node->next->key );
                            free( current_node->next );
                        }
                    free( table->table_data[ index ]->key );
                    free( table->table_data[ index ] );
                    table->table_data[ index ] = NULL;
                }
        }

    table->size = 0;
    free( table->table_data );
}


uint32_t generate_hash( const void *key,  int len, uint32_t seed )
{
 const unsigned int m = 0xc6a4a793;

  const int r = 16;

  unsigned int h = seed ^ (len * m);

  //----------
  
  const unsigned char * data = (const unsigned char *)key;

  while(len >= 4)
  {
    unsigned int k = *(unsigned int *)data;

    h += k;
    h *= m;
    h ^= h >> 16;

    data += 4;
    len -= 4;
  }
  
  //----------
  
  switch(len)
  {
  case 3:
    h += data[2] << 16;
    break;
  case 2:
    h += data[1] << 8;
    break;
  case 1:
    h += data[0];
    h *= m;
    h ^= h >> r;
    break;
  };
 
  //----------

  h *= m;
  h ^= h >> 10;
  h *= m;
  h ^= h >> 17;

  return h;

}

int ht_add( hash_table_t* table, char* to_add, void* add_val )
{
    uint32_t item_index;

    HT_Entry *new_entry = malloc( sizeof( HT_Entry ) );
    HT_Entry *current_node;

    new_entry->key = malloc( strlen( to_add ) + 1 );
    strcpy( new_entry->key, to_add );
    new_entry->value = add_val;

    new_entry->next = NULL;
    new_entry->prev = NULL;

    item_index = generate_hash( to_add, strlen( to_add ), HASH_NUMBER ) % table->capacity;

    // item not already in table
    if( table->table_data[ item_index ] == NULL )
        {
            table->table_data[ item_index ] = new_entry;
        }
    else
        {
            // item hash already in table
            current_node = table->table_data[ item_index ];
            while( current_node->next != NULL )
                {
                    // we don't want to add duplicates
                    if( strcmp( current_node->key, to_add ) != 0 )
                        {
                            // update the value
                            current_node = current_node->next;
                        }
                    else
                        {
                            free( new_entry->key );
                            free( new_entry );
                            return 0;
                        }
                }

            current_node->next = new_entry;
            new_entry->prev = current_node;

        }
    table->size++;

    return 1;
}


HT_Entry* find_item( hash_table_t* table, char* in_key )
{
    uint32_t search_index = generate_hash( in_key, strlen( in_key ), HASH_NUMBER ) % table->capacity;


    HT_Entry* current_node;

    if( ( table->size > 0 ) && ( table->table_data[ search_index ] != NULL ) )
        {

            current_node = table->table_data[ search_index ];

            while( strcmp( current_node->key, in_key ) != 0 ) 
                {

                    if( current_node->next == NULL )
                        {
                            return NULL;
                        }
                    current_node = current_node->next;
                }

            return current_node;
        }

    return NULL;
}


void *ht_find( hash_table_t* table, char* in_key )
{
    HT_Entry* found_item = find_item( table, in_key );

    if( found_item != NULL )
        {
            return found_item->value;
        }
    return NULL;
}

void* ht_delete( hash_table_t* table, char* in_key )
{
    HT_Entry *found_node = find_item( table, in_key );
    void* return_val;
    uint32_t found_index; 

    if( found_node != NULL && found_node->key != NULL )
        {
            return_val = found_node->value;
            found_index = generate_hash( in_key, strlen( in_key ), HASH_NUMBER ) % table->capacity;
            if( table->table_data[ found_index ] == found_node )
                {
                    table->table_data[ found_index ] = found_node->next;
                }
            if( found_node->next != NULL )
                {
                    found_node->next->prev = found_node->prev;
                }
            if( found_node->prev != NULL )
                {
                    found_node->prev->next = found_node->next;
                }

            free( found_node );

            table->size -= 1;

            return return_val;
        }


    return NULL;
}


HT_Entry *ht_get_items( hash_table_t* input )
{
    HT_Entry *output = NULL;
    HT_Entry* next_node;

    uint32_t input_index;
    uint32_t output_index;

    uint32_t capacity = input->capacity;

    if( input->size > 0 )
        {
            output = malloc( sizeof( HT_Entry ) * input->size );
            output_index = 0;

            for( input_index = 0; input_index < capacity; input_index++ )
                {
                    if( input->table_data[ input_index ] )
                        {
                            output[ output_index ] = *(input->table_data[ input_index ]);
                            next_node = input->table_data[ input_index ]->next;

                            output_index++;

                            while( next_node )
                                {
                                    output[ output_index ] = *(next_node);
                                    next_node = next_node->next;

                                    output_index++;
                                }
                        }
                }
        }
    return output;
}

HT_Entry *ht_get_items_no_malloc( hash_table_t* input, HT_Entry *data )
{
    HT_Entry *output = NULL;
    HT_Entry* next_node;

    uint32_t input_index;
    uint32_t output_index;

    uint32_t capacity = input->capacity;

    if( input->size > 0 )
        {
            output_index = 0;

            for( input_index = 0; input_index < capacity; input_index++ )
                {
                    if( input->table_data[ input_index ] )
                        {
                            data[ output_index ] = *(input->table_data[ input_index ]);
                            next_node = input->table_data[ input_index ]->next;

                            output_index++;

                            while( next_node )
                                {
                                    data[ output_index ] = *(next_node);
                                    next_node = next_node->next;

                                    output_index++;
                                }
                        }
                }
        }
    return output;
}
