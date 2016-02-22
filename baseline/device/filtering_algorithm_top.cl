/**********************************************************************
 * Nadesh Ramanathan, Imperial College London
 *
 * File: filtering_algorithm_top.cl 
 *
 * Additional Comments: distributed under a BSD license, see LICENSE.txt
 *
 **********************************************************************/

#include "filtering_algorithm_top.h"

typedef struct {
    node_pointer u;
    centre_list_pointer c;
    centre_index_type k;
    bool d;
} stack_record_type;

void init_allocator(uint par, local centre_list_pointer *flist, centre_list_pointer* next_free_location, const centre_list_pointer heapsize)
{
    for (centre_list_pointer i=0; i<=heapsize; i++) {
        flist[(uint)i] = i+1;
        if (i==heapsize) {
            break;
        }
    }
    *next_free_location = 1;
}

void k_malloc(uint par, local centre_list_pointer * flist, centre_list_pointer* next_free_location, centre_list_pointer *result)
{
    centre_list_pointer address = *next_free_location;
    *next_free_location = flist[(uint)address];
    *result = address;
}

void k_free(uint par, local centre_list_pointer * flist, centre_list_pointer* next_free_location, centre_list_pointer address)
{
    flist[(uint)address] = *next_free_location;
    *next_free_location = address;
}

void c_set_coord_type_vector_item(coord_type_vector *a, const coord_type b, const uint idx)
{
    a->value[idx] = b;
}


void c_set_coord_type_vector_ext_item(coord_type_vector_ext *a, const coord_type_ext b, const uint idx)
{
    a->value[idx] = b;
}

void c_set_coord_type_vector_ext_item_local(local coord_type_vector_ext *a, const coord_type_ext b, const uint idx)
{
    a->value[idx] = b;
}

void c_get_coord_type_vector_item(const coord_type_vector a, const uint idx, coord_type *result)
{
    coord_type tmp= a.value[idx];
    *result = tmp;
}


void c_get_coord_type_vector_ext_item(const coord_type_vector_ext a, const uint idx, coord_type_ext *result)
{
    coord_type_ext tmp= a.value[idx];
    *result = tmp;
}


void c_set_kd_tree_type_items(kdTree_type *a,
        const uint idx,
        const data_type_ext wgtCent,
        const data_type midPoint,
        const data_type bnd_hi,
        const data_type bnd_lo,
        const coord_type_ext sum_sq,
        const coord_type_ext count,
        const node_pointer left,
        const node_pointer right)
{
    a->idx        = idx;
    a->wgtCent    = wgtCent;
    a->midPoint   = midPoint; 
    a->bnd_hi     = bnd_hi;
    a->bnd_lo     = bnd_lo;
    a->sum_sq     = sum_sq;
    a->count      = count;
    a->left       = left;
    a->right      = right;
}


void c_get_kd_tree_type_items(local kdTree_type * restrict a,
        uint *idx,
        data_type_ext *wgtCent,
        data_type *midPoint,
        data_type *bnd_hi,
        data_type *bnd_lo,
        coord_type_ext *sum_sq,
        coord_type_ext *count,
        node_pointer *left,
        node_pointer *right
        )
{

    *idx        = a->idx;
    *wgtCent 	= a->wgtCent;
    *midPoint 	= a->midPoint;
    *bnd_hi		= a->bnd_hi;
    *bnd_lo		= a->bnd_lo;
    *sum_sq 	= a->sum_sq;
    *count 		= a->count;
    *left 		= a->left;
    *right 		= a->right;

}

// conversion from data_type to data_type_ext
void c_conv_short_to_long(data_type p, data_type_ext *result_out)
{
    data_type_ext result;
    for (uint d=0; d<D; d++) {
        coord_type tmp_r;
        c_get_coord_type_vector_item(p.value,d,&tmp_r);
        coord_type_ext tmp = (coord_type_ext) tmp_r; 
        c_set_coord_type_vector_ext_item(&result.value,tmp,d);
    }
    *result_out = result;
}

void c_saturate_mul_input(coord_type_ext val, mul_input_type *result)
{
    *result = (mul_input_type)val;
}


// fixed-point multiplication with saturation and scaling
void c_fi_mul(coord_type_ext op1, coord_type_ext op2, coord_type_ext *result)
{
    mul_input_type tmp_op1;
    c_saturate_mul_input(op1,&tmp_op1);
    mul_input_type tmp_op2;
    c_saturate_mul_input(op2,&tmp_op2);

    long int result_unscaled;
    result_unscaled = tmp_op1*tmp_op2;

    long int result_scaled;
    result_scaled = result_unscaled >> MUL_FRACTIONAL_BITS;
    *result = (coord_type_ext)result_scaled;
}


// tree adder
void c_tree_adder(coord_type_ext *input_array, coord_type_ext *result)
{
    coord_type_ext tmp = 0;
    for(uint d=0; d<D; d++){
        tmp += input_array[d];
    }
    *result = tmp;
}

// tree adder (overloaded function)
void c_tree_adder_m(local coord_type_ext *input_array,const uint m, coord_type_ext *result)
{
    coord_type_ext tmp = 0;
    for(uint p=0; p<P; p++){
        tmp += input_array[p];
    }
    *result =tmp;

}

// inner product of p1 and p2
void inline c_dot_product(data_type_ext p1,data_type_ext p2, coord_type_ext *r)
{

    coord_type_ext tmp_mul_res[D];

    for (uint d=0;d<D;d++) {
        mul_input_type tmp_op1;
        coord_type_ext c_tmp1;
        c_get_coord_type_vector_ext_item(p1.value,d,&c_tmp1);
        c_saturate_mul_input(c_tmp1,&tmp_op1);

        mul_input_type tmp_op2;
        coord_type_ext c_tmp2;
        c_get_coord_type_vector_ext_item(p2.value,d,&c_tmp2);
        c_saturate_mul_input(c_tmp2,&tmp_op2);
        coord_type_ext tmp_mul = tmp_op1*tmp_op2;
        tmp_mul_res[d] = tmp_mul;
    }

    coord_type_ext tmp;
    c_tree_adder(tmp_mul_res,&tmp);
    *r = tmp;
}

// compute the Euclidean distance between p1 and p2
void inline c_compute_distance(data_type_ext p1, data_type_ext p2, coord_type_ext *dist)
{

    data_type_ext tmp_p1 = p1;
    data_type_ext tmp_p2 = p2;
    coord_type_ext tmp_mul_res[D];

    for (uint d=0; d<D; d++) {
        coord_type_ext tmp_sub1;
        c_get_coord_type_vector_ext_item(tmp_p1.value,d,&tmp_sub1);
        coord_type_ext tmp_sub2;
        c_get_coord_type_vector_ext_item(tmp_p2.value,d,&tmp_sub2);
        coord_type_ext tmp = tmp_sub1 - tmp_sub2;
        coord_type_ext tmp_mul; 
        c_fi_mul(tmp,tmp,&tmp_mul);
        tmp_mul_res[d] = tmp_mul;
    }

    coord_type_ext tmp;
    c_tree_adder(tmp_mul_res,&tmp);
    *dist = tmp;
}

// check whether any point of bounding box is closer to z than to z*
// this is a modified version of David Mount's code (http://www.cs.umd.edu/~mount/Projects/KMeans/ )
void inline c_tooFar_fi(data_type closest_cand, data_type cand, data_type bnd_lo, data_type bnd_hi, bool *too_far)
{

    coord_type_ext boxDot;
    coord_type_ext ccDot;

    data_type_ext tmp_closest_cand;
    c_conv_short_to_long(closest_cand,&tmp_closest_cand);

    data_type_ext tmp_cand;
    c_conv_short_to_long(cand,&tmp_cand);

    data_type_ext tmp_bnd_lo;
    c_conv_short_to_long(bnd_lo,&tmp_bnd_lo);

    data_type_ext tmp_bnd_hi;
    c_conv_short_to_long(bnd_hi,&tmp_bnd_hi);

    coord_type_ext tmp_mul_res[D];
    coord_type_ext tmp_mul_res2[D];

    for (uint d = 0; d<D; d++) {
        coord_type_ext tmp_sub_op1; 
        c_get_coord_type_vector_ext_item(tmp_cand.value,d,&tmp_sub_op1);
        coord_type_ext tmp_sub_op2;
        c_get_coord_type_vector_ext_item(tmp_closest_cand.value,d,&tmp_sub_op2);
        coord_type_ext ccComp = tmp_sub_op1-tmp_sub_op2;
        coord_type_ext tmp_mul;
        c_fi_mul(ccComp,ccComp,&tmp_mul);
        tmp_mul_res[d] = tmp_mul;

        coord_type_ext tmp_diff2;
        coord_type_ext tmp_sub_op3;
        coord_type_ext tmp_sub_op4;
        if (ccComp > 0) {
            c_get_coord_type_vector_ext_item(tmp_bnd_hi.value,d,&tmp_sub_op3);
        }
        else {
            c_get_coord_type_vector_ext_item(tmp_bnd_lo.value,d,&tmp_sub_op3);
        }
        c_get_coord_type_vector_ext_item(tmp_closest_cand.value,d,&tmp_sub_op4);
        tmp_diff2 = tmp_sub_op3 - tmp_sub_op4;

        coord_type_ext tmp_mul2;
        c_fi_mul(tmp_diff2,ccComp,&tmp_mul2);
        tmp_mul_res2[d] = tmp_mul2;

    }

    c_tree_adder(tmp_mul_res,&ccDot);
    c_tree_adder(tmp_mul_res2,&boxDot);

    coord_type_ext tmp_boxDot = boxDot<<1;
    bool tmp_res;
    if (ccDot>tmp_boxDot) {
        tmp_res = true;
    } else {
        tmp_res = false;
    }

    *too_far = tmp_res;

}

void c_init_stack(uint *stack_pointer)
{
    *stack_pointer = 0;

}

// push pointer to tree node pointer onto stack
void c_push(uint par,node_pointer u, centre_list_pointer c, centre_index_type k,  bool d, uint *stack_pointer, local stack_record_type *stack_array, uint *result )
{
    uint tmp = *stack_pointer;

    stack_array[tmp].u = u;
    stack_array[tmp].c = c;
    stack_array[tmp].k = k;
    stack_array[tmp].d = d;

    tmp++;
    *stack_pointer=tmp;
    *result = tmp; 
    // return tmp;
}

// pop pointer to tree node pointer from stack
void c_pop(uint par,node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, local stack_record_type *stack_array, uint *result)
{
    uint tmp = *stack_pointer-1;
    *u = stack_array[tmp].u;
    *c = stack_array[tmp].c;
    *k = stack_array[tmp].k;
    *d = stack_array[tmp].d;
    *stack_pointer = tmp;
    *result = tmp;
    //    return tmp;
}

// look up head of node stack
void c_lookahead(uint par, node_pointer *u, centre_list_pointer *c, centre_index_type *k, bool *d, uint *stack_pointer, local stack_record_type *stack_array, uint *result )
{
    uint tmp = *stack_pointer-1;
    *u = stack_array[tmp].u;
    *c = stack_array[tmp].c;
    *k = stack_array[tmp].k;
    *d = stack_array[tmp].d;
    *result = tmp;
    //    return tmp;
}


// update the new centre positions after one outer clustering iteration
void update_centres(local centre_type * restrict centres_in,centre_index_type k, uint batchsize, local data_type * restrict centres_positions_out)
{
    uint lid = get_local_id(0);
    for (centre_index_type i=0; i<(batchsize); i++) {
        uint addr = lid * (batchsize) + i;
        coord_type_ext tmp_count = centres_in[addr].count;
        if ( tmp_count == 0 )
            tmp_count = 1;

        data_type_ext tmp_wgtCent = centres_in[addr].wgtCent;
        data_type tmp_new_pos;
        for (uint d=0; d<D; d++) {
            coord_type_ext d_tmp;
            c_get_coord_type_vector_ext_item(tmp_wgtCent.value,d,&d_tmp);
            coord_type_ext tmp_div_ext = ( d_tmp / tmp_count); //let's see what it does with that.
            coord_type tmp_div = (coord_type) tmp_div_ext;
            c_set_coord_type_vector_item(&tmp_new_pos.value,tmp_div,d);
        }
        centres_positions_out[addr] = tmp_new_pos;
    }
}


void prune_centre_set( 
        local centre_heap_type *restrict heap_in,
        global centre_heap_type *restrict heap_out,
        centre_index_type k,
        //centre_index_type *centre_indices_in,
        local data_type *centre_positions,
        data_type u_bnd_hi,
        data_type u_bnd_lo,
        data_type z_star,
        //centre_index_type *centre_set_out,
        centre_index_type *k_out
        )
{

    //copy candidates that survive pruning into new list
    centre_index_type new_k=0; 
    centre_index_type tmp_new_idx=0;


    // determine whether a sub-tree will be pruned
    for (centre_index_type i=0; i<=k; i++) {
        bool too_far;

        centre_index_type tmp_index = heap_in->idx[i];

        data_type position = centre_positions[tmp_index];
        c_tooFar_fi(z_star, position, u_bnd_lo, u_bnd_hi, &too_far);
        if ( too_far==false ) {

            heap_out->idx[tmp_new_idx] = tmp_index;

#ifdef VERBOSE
#ifndef SYN
            printf("%d ", tmp_index);
#endif
#endif

            tmp_new_idx++;
            new_k++;
        }
        if (i==k) {
            break;
        }
    }

    *k_out = new_k-1;
}

void minsearch( 
        local centre_heap_type *restrict heap_in,
        centre_index_type k,
        //centre_index_type *centre_set,
        data_type_ext comp_point,
        local data_type *centre_positions,
        centre_index_type *min_index,
        data_type *z_star
        //centre_index_type *centre_indices_out
        )
{

    centre_index_type tmp_final_idx;
    data_type tmp_z_star;
    coord_type_ext tmp_min_dist;
    tmp_min_dist = (1<<(COORD_BITWITDH_EXT-1))-1;

    // find centre with smallest distance to z_star
    for (centre_index_type i=0; i<=k; i++) {

        centre_index_type tmp_index = heap_in->idx[i];

#ifdef VERBOSE
#ifndef SYN
        printf("%d ", tmp_index);
#endif
#endif

        coord_type_ext tmp_dist;
        data_type position = centre_positions[tmp_index];

        data_type_ext tmp_position;
        c_conv_short_to_long(position,&tmp_position);
        c_compute_distance(tmp_position, comp_point, &tmp_dist);

        if ((tmp_dist < tmp_min_dist) ) {
            tmp_min_dist = tmp_dist;
            tmp_final_idx = tmp_index;
            tmp_z_star = position;
        }

    }

#ifdef VERBOSE
#ifndef SYN
    printf(" # ");
#endif
#endif

    *z_star = tmp_z_star;
    *min_index = tmp_final_idx;
}

void calculate_distortion(data_type_ext u_wgtCent,
        coord_type_ext u_count,
        coord_type_ext u_sum_sq,
        data_type z_star,
        coord_type_ext *sum_sq_out)
{
    // some scaling...
    data_type_ext tmp_wgtCent = u_wgtCent;
    for (uint d=0; d<D; d++) {
        coord_type_ext tmp;
        c_get_coord_type_vector_ext_item(tmp_wgtCent.value,d,&tmp);
        tmp = tmp >> MUL_FRACTIONAL_BITS;
        c_set_coord_type_vector_ext_item(&tmp_wgtCent.value,tmp,d);
    }

    // update sum_sq of centre
    coord_type_ext tmp1_2, tmp2_2;
    data_type_ext tmp_z_star;
    c_conv_short_to_long(z_star,&tmp_z_star);

    c_dot_product(tmp_z_star,tmp_wgtCent,&tmp1_2);
    c_dot_product(tmp_z_star,tmp_z_star ,&tmp2_2);
    coord_type_ext tmp1, tmp2;
    tmp1 = tmp1_2<<1;
    tmp2 = tmp2_2>>1;
    coord_type_ext tmp_count = u_count;
    coord_type_ext tmp2_sat;
    c_saturate_mul_input(tmp2,&tmp2_sat);
    coord_type_ext tmp_count_sat;
    c_saturate_mul_input(tmp_count,&tmp_count_sat);
    coord_type_ext tmp3 = tmp2_sat*tmp_count_sat;
    coord_type_ext tmp_sum_sq1 = u_sum_sq+tmp3;
    coord_type_ext tmp_sum_sq = tmp_sum_sq1-tmp1;
    *sum_sq_out = tmp_sum_sq;
}

void process_node( uint par,
        local centre_heap_type * restrict heap_in,
        global centre_heap_type * restrict heap_out,
        centre_index_type k,
        data_type_ext u_wgtCent,
        data_type u_midPoint,
        data_type u_bnd_hi,
        data_type u_bnd_lo,
        coord_type_ext u_sum_sq,
        coord_type_ext u_count,
        node_pointer u_left,
        node_pointer u_right,
        //centre_index_type *centre_set_data,
        local data_type *centre_positions,
        centre_index_type *k_out,
        //centre_index_type *centre_indices_out,
        centre_index_type *final_centre_index,
        coord_type_ext *sum_sq_out,
        uint *dead_end )
{
    centre_index_type tmp_k = k;

    uint tmp_deadend;
    centre_index_type tmp_final_centre_index;
    coord_type_ext tmp_sum_sq_out;
    centre_index_type tmp_k_out;

    // leaf node?
    data_type_ext comp_point;
    if ( (u_left == NULL_PTR) && (u_right == NULL_PTR) ) {
        comp_point = u_wgtCent;
    } else {
        data_type_ext tmp_midPoint;
        c_conv_short_to_long(u_midPoint,&tmp_midPoint);
        comp_point = tmp_midPoint;
    }

    centre_index_type tmp_final_idx;
    data_type z_star;

    minsearch( 
            heap_in,
            tmp_k,
            //centre_set_data,
            comp_point,
            centre_positions,
            &tmp_final_idx,
            &z_star
            //tmp_centre_indices
            );

    centre_index_type new_k;
    coord_type_ext tmp_sum_sq;

    prune_centre_set( 
            heap_in,
            heap_out,
            tmp_k,
            //tmp_centre_indices,
            centre_positions,
            u_bnd_hi,
            u_bnd_lo,
            z_star,
            //centre_indices_out,
            &new_k );

    calculate_distortion( u_wgtCent,
            u_count,
            u_sum_sq,
            z_star,
            &tmp_sum_sq);

    if ((new_k == 0) || ( (u_left == NULL_PTR) && (u_right == NULL_PTR) )) {
        tmp_deadend = 1;
    } else {
        tmp_deadend = 0;
    }

    *k_out = new_k;
    *final_centre_index = tmp_final_idx;
    *sum_sq_out = tmp_sum_sq;
    *dead_end = tmp_deadend;

}


// main clustering kernel
void filter (
        global kdTree_type (* restrict tree)[TREE_HEAP_SIZE],
        local kdTree_type * restrict tree_node,
        global node_pointer * restrict root_array,
        const centre_index_type k,
        local centre_list_pointer * restrict centre_freelist,
        local data_type * restrict centre_positions,
        local centre_type * restrict centre_buffer,
        local stack_record_type * restrict stack_array,
        local centre_heap_type * restrict private_heap,
        global centre_heap_type * restrict cntr_heap,
        global uint * restrict vnodes,
        global uint * restrict ncpairs,
        global uint * restrict allocs
        )
{

    uint par = get_local_id(0);
    uint visited_nodes = 0;
    uint nc_pairs = 0;
    uint allocated_centre_sets = 1;
    uint max_allocated_centre_sets = 1;

    node_pointer root = root_array[par];

#ifdef CENTRE_BUFFER_ONCHIP
    // init centre buffer
    for(centre_index_type i=0; i<k ; i++) {
        centre_type tmp;
        tmp.count = 0;
        tmp.sum_sq = 0;
        tmp.wgtCent.value.value[0] = 0;
        tmp.wgtCent.value.value[1] = 0;
        tmp.wgtCent.value.value[2] = 0;

        centre_buffer[i] = tmp;
    }
#endif

    // dynamic memory allocation
    centre_list_pointer centre_next_free_location;

    // init dynamic memory allocator for centre lists scratchpad heap
    init_allocator(par,centre_freelist, &centre_next_free_location, CENTRESET_HEAP_SIZE-2);

    // stack pointer
    uint stack_pointer;
    c_init_stack(&stack_pointer);
    uint node_stack_length;

    centre_list_pointer centre_list_idx;
    k_malloc(par,centre_freelist, &centre_next_free_location, &centre_list_idx);

    // write the malloc'ed data set to the first valid address
    for(centre_index_type i=0; i<=k; i++) {
        cntr_heap[centre_list_idx].idx[i] = i;
    }


    // push pointers to tree root node and first centre list onto the stack
    c_push(par,root, centre_list_idx, k,  true, &stack_pointer, stack_array, &node_stack_length);

    // main tree search loop
    while (node_stack_length != 0) {


        // fetch head of stack
        node_pointer u;
        centre_list_pointer centre_set_in,centre_set_out;
        centre_index_type tmp_k;
        bool rdy_for_deletion;



        // fetch pointer to centre list
        c_pop(par,&u, &centre_set_in, &tmp_k, &rdy_for_deletion, &stack_pointer, stack_array, &node_stack_length);

        // fetch tree node
        *tree_node = tree[par][(uint)(u)];


        uint dummy_idx;
        data_type_ext u_wgtCent;
        data_type u_midPoint;
        data_type u_bnd_lo;
        data_type u_bnd_hi;
        coord_type_ext u_sum_sq;
        coord_type_ext u_count;
        node_pointer u_left;
        node_pointer u_right;
        c_get_kd_tree_type_items(	tree_node,
                &dummy_idx,
                &u_wgtCent,
                &u_midPoint,
                &u_bnd_hi,
                &u_bnd_lo,
                &u_sum_sq,
                &u_count,
                &u_left,
                &u_right);


#ifdef VERBOSE
#ifndef SYN
        printf("(%d, %d, %d); ",u_count, u_left, u_right);
#endif

#ifndef SYN
        printf("R:[%d] ", centre_set_in);
#endif
#endif


        k_malloc(par,centre_freelist, &centre_next_free_location,&centre_set_out);
#ifndef SYN
        allocated_centre_sets++;
        if (max_allocated_centre_sets < allocated_centre_sets)
            max_allocated_centre_sets = allocated_centre_sets;
#endif

        centre_index_type tmp_k_out;
        centre_heap_type tmp_centre_indices_out;
        centre_index_type tmp_final_centre_index;
        coord_type_ext tmp_sum_sq_out;
        uint tmp_deadend;

        for (centre_index_type i=0; i<=tmp_k; i++) {
            private_heap->idx[i] = cntr_heap[centre_set_in].idx[i];
        }

        process_node(       par,
                private_heap,    
                &cntr_heap[centre_set_out],
                tmp_k,
                u_wgtCent,
                u_midPoint,
                u_bnd_hi,
                u_bnd_lo,
                u_sum_sq,
                u_count,
                u_left,
                u_right,
                //tmp_cntr_list_data.idx,
                centre_positions,
                &tmp_k_out,
                //tmp_centre_indices_out.idx,
                &tmp_final_centre_index,
                &tmp_sum_sq_out,
                &tmp_deadend );

        // free list that has been read twice
        if (rdy_for_deletion == true) {
            k_free(par,centre_freelist, &centre_next_free_location, centre_set_in);
#ifndef SYN
            allocated_centre_sets--;
#endif
        }

        // write back
        // final decision whether sub-tree will be pruned
        if ( tmp_deadend == 1 ) {

            k_free(par,centre_freelist, &centre_next_free_location, centre_set_out);
#ifndef SYN
            allocated_centre_sets--;
#endif
update_function: {
                     // weighted centroid of this centre
                     centre_type tmpCentreRecord = centre_buffer[tmp_final_centre_index];
                     for (uint d=0; d<D; d++) {
                         coord_type_ext tmp1;
                         c_get_coord_type_vector_ext_item(tmpCentreRecord.wgtCent.value,d,&tmp1);
                         coord_type_ext tmp2;
                         c_get_coord_type_vector_ext_item(u_wgtCent.value,d,&tmp2);
                         c_set_coord_type_vector_ext_item(&tmpCentreRecord.wgtCent.value,tmp1+tmp2,d);
                     }
                     // update number of points assigned to centre
                     coord_type_ext tmp1 =  u_count;
                     coord_type_ext tmp2 =  tmpCentreRecord.count;
                     tmpCentreRecord.count = tmp1 + tmp2;
                     coord_type_ext tmp3 =  tmp_sum_sq_out;
                     coord_type_ext tmp4 =  tmpCentreRecord.sum_sq;
                     tmpCentreRecord.sum_sq  = tmp3 + tmp4;
                     centre_buffer[tmp_final_centre_index] = tmpCentreRecord;
                 }

        } else {

            // allocate new centre list
            centre_index_type new_k = tmp_k_out;

            node_pointer left_child = u_left;
            node_pointer right_child = u_right;

#ifdef VERBOSE
#ifndef SYN
            printf("W:[%d] ", centre_set_out);
#endif
#endif

            // push children onto stack
            c_push(par,right_child, centre_set_out,new_k, true, &stack_pointer, stack_array, &node_stack_length);
            c_push(par,left_child, centre_set_out,new_k, false, &stack_pointer, stack_array, &node_stack_length);

        }

#ifdef VERBOSE
#ifndef SYN
        printf("\n");
#endif
#endif

        visited_nodes++;
        nc_pairs += tmp_k;

    }

#ifndef SYN
    printf("%d: node-center pairs: %d\n",par,nc_pairs);
    printf("%d: visited nodes: %d\n",par,visited_nodes);
    printf("%d: max allocated centre sets: %d\n",par,max_allocated_centre_sets);
#endif

    vnodes[par] = visited_nodes;
    ncpairs[par] = nc_pairs;
    allocs[par] = allocated_centre_sets;


}

__attribute__((reqd_work_group_size(P,1,1)))
    kernel void filtering_algorithm_top ( 
            global kdTree_type (* restrict tree_image)[TREE_HEAP_SIZE],   
            global data_type * restrict cntr_pos_init, 
            global node_pointer * restrict root, 
            constant centre_index_type * restrict k_in,
            global centre_heap_type (*restrict cntr_heap)[CENTRESET_HEAP_SIZE],
            global coord_type_ext * restrict distortion_out,
            global data_type * restrict clusters_out,
            global uint (*restrict vnodes_out)[P],
            global uint (*restrict ncpairs_out)[P],
            global uint (*restrict alloc_out)[P]
            ){

        const uint gid = get_global_id(0);
        const uint lid = get_local_id(0);
        const uint par = get_local_size(0);
        const uint k = *k_in;
        const uint k_batchsize = k/P;

        //localized inputs 
        local centre_type filt_centres_out_reduced[K];
        local data_type new_centre_positions[K];

        // stack, heap, tree node, heap list and result buffers
        local stack_record_type stack_array[P][STACK_SIZE];
        local centre_heap_type private_heap[P];
        local kdTree_type tree_node[P];
        local centre_list_pointer centre_freelist[P][CENTRESET_HEAP_SIZE];
        local centre_type centre_buffer[P][K];

        // Reduction variables
        local coord_type_ext arr_count[K][P];
        local coord_type_ext arr_sum_sq[K][P];
        local coord_type_vector_ext arr_wgtCent[K][P];
        local coord_type_vector_ext tmp_sum[K];
        local coord_type_ext tmp_a[K][D][P];
        barrier(CLK_LOCAL_MEM_FENCE);


        if(lid == 0){
            for (uint i=0; i<k ; i++) {
                data_type position;
                position = cntr_pos_init[i];
                new_centre_positions[i] = position;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);


        // iterate over a constant number of outer clustering iterations
        for (uint l=0; l<L; l++) {

            barrier(CLK_LOCAL_MEM_FENCE);

            filter(tree_image, &tree_node[lid], root, k-1, centre_freelist[lid], new_centre_positions, centre_buffer[lid], stack_array[lid],&private_heap[lid],cntr_heap[lid],vnodes_out[l],ncpairs_out[l],alloc_out[l]); 

            barrier(CLK_LOCAL_MEM_FENCE);

#if P>1
            for(centre_index_type i=0; i < k; i++) {
                arr_count[i][lid] = ((coord_type_ext)centre_buffer[lid][i].count);
                arr_sum_sq[i][lid] = (centre_buffer[lid][i].sum_sq);
                arr_wgtCent[i][lid] = (centre_buffer[lid][i].wgtCent.value);

                for (uint d=0; d<D; d++) {
                    coord_type_ext tmp_h;
                    tmp_h = arr_wgtCent[i][lid].value[d];
                    tmp_a[i][d][lid] = tmp_h;
                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);
            for(centre_index_type i=0; i < k_batchsize; i++) {
                uint addr = lid*(k_batchsize) + i;
                coord_type_ext tmp1 = 0;
                coord_type_ext tmp2 = 0;

                c_tree_adder_m(arr_count[addr],P,&tmp1);
                c_tree_adder_m(arr_sum_sq[addr],P,&tmp2);
                filt_centres_out_reduced[addr].count = tmp1;
                filt_centres_out_reduced[addr].sum_sq = tmp2;

                for (uint d=0; d<D; d++) {
                    coord_type_ext tmp = 0;
                    c_tree_adder_m(tmp_a[addr][d],P,&tmp);
                    c_set_coord_type_vector_ext_item_local(&tmp_sum[addr],tmp,d);
                }
                filt_centres_out_reduced[addr].wgtCent.value = tmp_sum[addr];
            }
#else

            for(centre_index_type i=0; i < k; i++){
                filt_centres_out_reduced[i] = centre_buffer[0][i];
            }
#endif

            barrier(CLK_LOCAL_MEM_FENCE);

            // re-init centre positions
            update_centres(filt_centres_out_reduced, k, k_batchsize, new_centre_positions);

            barrier(CLK_LOCAL_MEM_FENCE);

#ifndef SYN
            if(lid==0){
                printf("Round %d Centreset Results\n",l);
                for (uint i=0; i<k; i++) {
                    printf("%d: ",i);
                    for (uint d=0; d<D-1; d++) {
                        printf("%d ",new_centre_positions[i].value.value[d]);
                    }
                    printf("%d, ",new_centre_positions[i].value.value[D-1]);

                    printf("%d\n",filt_centres_out_reduced[i].sum_sq);

                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);

#endif
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        for (centre_index_type addr=0; addr<(k_batchsize); addr++) {
            uint i = lid*(k_batchsize) + addr;
            distortion_out[i] = filt_centres_out_reduced[i].sum_sq;
            clusters_out[i]   = new_centre_positions[i];
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
