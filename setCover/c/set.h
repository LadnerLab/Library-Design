#ifndef SET_H_INCLUDED
#define SET_H_INCLUDED
#include "hash_table.h"

typedef struct set_t
{
    hash_table_t* data;
} set_t;

/**
 * Initializes a set
 * Initializes set's hash_table data
 * @param init set member to initialize
 **/
void set_init( set_t* to_init );

/**
 * Adds a string value to a set.
 * Note: Adds item to set's hash table
 * @param set_to_add Set to add member to
 * @param add_data character string to add to
 **/
void set_add( set_t* set_to_add, char* add_data );
int set_check( set_t* source, char* item );

/**
 * Removes a string from a set
 * Note: Removes object from internal hash_table,
 *       which frees the key stored there
 * @param set_to_remove set to remove string data from
 * @param remove_data string data to remove
 * @returns integer boolean success of operation 
 **/
int set_remove( set_t* set_to_remove, char* remove_data );

/**
 * clears a set object and frees all of its memory
 * and the memory of its members
 * @param set_to_clear set to remove data from
 **/
void set_clear( set_t* set_to_clear );

/**
 * Adds all of the elements found in source to dest
 * @param dest set destination to add items to
 * @param source set to add items to dest
 **/ 
void set_update( set_t* dest, set_t* source );


/**
 * Computes the relative complement of sets first and second
 * e.g. x an element of first where x is not an element of second
 * Difference is stored in dest set
 * @param dest set to store results in
 * @param first set to compare difference of
 * @param second set to compare difference of
 **/
void set_difference( set_t* first, set_t* second );



/**
 * Adds num_elements items from an array of strings into a set
 * @param dest set to add items to
 * @param in_array string array of items to add
 * @param num_elements number of elements to add to set 
 **/
void set_add_all( set_t* dest, char** in_array, int num_elements );
#endif
