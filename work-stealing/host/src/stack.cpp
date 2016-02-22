/**********************************************************************
 * Nadesh Ramanathan, Imperial College London
 *
 * File: stack.cpp
 *
 * Additional Comments: distributed under a BSD license, see LICENSE.txt
 *
 **********************************************************************/


#include "stack.h"

stack_record_type& stack_record_type::operator=(const stack_record_type& a)
{
    u = a.u;
    c = a.c;
    k = a.k;
    d = a.d;
    return *this;
}


void init_stack(uint *stack_pointer)
{
    *stack_pointer = 0;

}


// push pointer to tree node pointer onto stack
uint push(node_pointer u, centre_list_pointer c, centre_index_type k,  bool d, uint *stack_pointer, stack_record_type* stack_array)
{
    uint tmp = *stack_pointer;

    stack_array[tmp].u = u;
    stack_array[tmp].c = c;
    stack_array[tmp].k = k;
    stack_array[tmp].d = d;

    tmp++;
    *stack_pointer=tmp;
    return tmp;
}

// pop pointer to tree node pointer from stack
uint pop(node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, stack_record_type* stack_array)
{
    uint tmp = *stack_pointer-1;
    *u = stack_array[tmp].u;
    *c = stack_array[tmp].c;
    *k = stack_array[tmp].k;
    *d = stack_array[tmp].d;
    *stack_pointer = tmp;
    return tmp;
}

// look up head of node stack
uint lookahead(node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, stack_record_type* stack_array)
{
    uint tmp = *stack_pointer-1;
    *u = stack_array[tmp].u;
    *c = stack_array[tmp].c;
    *k = stack_array[tmp].k;
    *d = stack_array[tmp].d;
    return tmp;
}


// push pointer to tree node pointer onto stack
uint push_on_chip(node_pointer u, centre_list_pointer c, centre_index_type k,  bool d, uint *stack_pointer, stack_record_type *stack_memory)
{
    uint tmp = *stack_pointer;

    stack_record_type in;
    in.u = u;
    in.c = c;
    in.k = k;
    in.d = d;

    stack_memory[tmp] = in;

    tmp++;
    *stack_pointer=tmp;
    return tmp;
}


// pop pointer to tree node pointer from stack
uint pop_on_chip(node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, stack_record_type *stack_memory)
{
    uint tmp = *stack_pointer-1;

    stack_record_type out;

    out = stack_memory[tmp];

    *u = out.u;
    *c = out.c;
    *k = out.k;
    *d = out.d;
    *stack_pointer = tmp;
    return tmp;
}


