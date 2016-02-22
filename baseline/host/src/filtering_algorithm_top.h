/**********************************************************************
 * Nadesh Ramanathan, Imperial College London
 *
 * File: filtering_algorithm_top.h
 *
 * Additional Comments: distributed under a BSD license, see LICENSE.txt
 *
 **********************************************************************/

#ifndef FILTERING_ALGORITHM_TOP_H
#define FILTERING_ALGORITHM_TOP_H

#define D 3         // data dimensionality
#define N 1024*1024 // max number of data points

#define P 4 
#define K 128       // max number of centres
#define L 16        // max number of outer clustering iterations

#define STACK_SIZE 128
#define TREE_HEAP_SIZE 4*N      // max size of heap memory for the kd-tree (2*n tree nodes)
#define CENTRESET_HEAP_SIZE 256 // max number of centre lists that can be allocated in the scratchpad heap

#define CENTRE_BUFFER_ONCHIP
#define COORD_BITWIDTH 16
#define COORD_BITWITDH_EXT 32
#define NODE_POINTER_BITWIDTH 32    // log2(2*N)

// pointer types to tree nodes and centre lists
typedef unsigned int uint;
typedef uint node_pointer;
typedef uint centre_list_pointer;
typedef uint centre_index_type;

typedef  short int coord_type;
typedef  int coord_type_ext;

typedef struct {coord_type value[D+1];} coord_type_vector ;
typedef struct { coord_type_ext value[D+1];} coord_type_vector_ext ;

// ... used for saturation
#define MAX_FIXED_POINT_VAL_EXT ((1<<(COORD_BITWITDH_EXT-1))-1)

//bit width definitions for multiplication
#define MUL_INTEGER_BITS 12
#define MUL_FRACTIONAL_BITS 6
#define MUL_MAX_VAL ((1<<(MUL_INTEGER_BITS+MUL_FRACTIONAL_BITS-1))-1)
#define MUL_MIN_VAL (-1*(1<<(MUL_INTEGER_BITS+MUL_FRACTIONAL_BITS-1)))
#define NULL_PTR 0

// this should be always 1
#define FILE_INDEX 1

#define TREE_NODE_BITWIDTH (3*D*COORD_BITWIDTH+D*COORD_BITWITDH_EXT+2*COORD_BITWITDH_EXT+2*NODE_POINTER_BITWIDTH+32)

typedef int mul_input_type;

typedef struct {
    coord_type_vector value;
} data_type;

typedef struct {
    coord_type_vector_ext value;
} data_type_ext;

typedef struct {
    uint            idx;
    data_type_ext   wgtCent;
    data_type       midPoint;
    data_type       bnd_hi;
    data_type       bnd_lo;
    coord_type_ext  sum_sq;
    coord_type_ext  count;
    node_pointer    left;
    node_pointer    right;
} kdTree_type ;

typedef kdTree_type* kdTree_ptr;

typedef struct {
    data_type_ext wgtCent; // sum of all points assigned to this centre
    coord_type_ext sum_sq; // sum of norm of all points assigned to this centre
    coord_type_ext count;
} centre_type;

typedef centre_type* centre_ptr;


typedef struct {
    centre_index_type idx[K];
}centre_heap_type;

#endif  /* FILTERING_ALGORITHM_TOP_H */
