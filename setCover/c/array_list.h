#ifndef ARRAY_LIST_H_INCLUDED
#define ARRAY_LIST_H_INCLUDED
#include <stdint.h>
typedef struct array_list_t
{
    uint32_t size;
    uint32_t capacity;

    void **array_data;
} array_list_t;

/**
 * Initializes an array_list_t object by setting size to zero,
 * and allocating capacity for array_data
 * @param to_init pointer to array_list_t object to initialize
 **/
void ar_init( array_list_t* to_init );

/**
 * clears an arraylist object, frees each object in 
 * array_data, and then frees to_clear itself
 * @param to_clear pointer to array_list_t object to clear
 **/
void ar_clear( array_list_t* to_clear );

/**
 * Clears an array list as specified by ar_clear,
 * but also frees each entry in its table 
 * 
 * @param to_clear pointer to array_list to clear
 **/
void ar_clear_and_free( array_list_t* to_clear );

/**
 * Retrieves item from arraylist at index
 * @param to_get pointer to array_list_t to retrieve item from
 * @param index unsigned integer index at which to retrieve
 * @returns void pointer to object found at index, or null otherwise
 **/
void *ar_get( array_list_t* to_get, unsigned int index );

/**
 * Adds an item to specified array_list
 * Note: Automatically resizes the array_list if necessary
 * @param to_add pointer to array_list_t object to add to
 * @new_data void data to add to array_list
 **/
void ar_add( array_list_t *to_add, void* new_data );

/**
 * Sets element in to_set at index to new data
 * Note: Index must be less than the size of the array_list
 * @param to_set to_set array_list pointer to set data at
 * @param index index of array_list's data to set 
 * @param new_data new data to add to array_list
 **/
void ar_set( array_list_t *to_set, unsigned int index, void* new_data );

/**
 * Remove an item from array list
 * Note: Index must be greater than zero and less than the size of
 *       the array list
 * @param to_remove pointer to array_list to remove an item from
 * @param remove_index index of item to remove
 * @returns pointer to item removed 
 **/
void *ar_remove( array_list_t *to_remove, unsigned int remove_index );


#endif
