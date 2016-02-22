/**********************************************************************
* Nadesh Ramanathan, Imperial College London
*
* File: stack.h
*
* Additional Comments: distributed under a BSD license, see LICENSE.txt
*
**********************************************************************/


#ifndef STACK_H
#define STACK_H

#include "filtering_algorithm_top.h"


// centre stack
struct stack_record_type {
	node_pointer u;
    centre_list_pointer c;
    centre_index_type k;
    bool d;
    stack_record_type& operator=(const stack_record_type& a);
};

#define WORDS_PER_STACK_RECORD 4

void init_stack(uint *stack_pointer);
uint push(node_pointer u, centre_list_pointer c, centre_index_type k,  bool d, uint *stack_pointer, stack_record_type* stack_array);
uint pop(node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, stack_record_type* stack_array);
uint lookahead(node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, stack_record_type* stack_array);


#endif
