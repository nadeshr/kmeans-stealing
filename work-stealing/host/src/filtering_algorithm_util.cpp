/**********************************************************************
* Nadesh Ramanathan, Imperial College London
*
* File: filtering_algorithm_util.cpp
*
* Additional Comments: distributed under a BSD license, see LICENSE.txt
*
**********************************************************************/

#include <math.h>
#include "filtering_algorithm_util.h"


/* ****** helper functions *******/

void set_coord_type_vector_item(coord_type_vector *a, const coord_type b, const uint idx)
{
    a->value[idx] = b;
}


void set_coord_type_vector_ext_item(coord_type_vector_ext *a, const coord_type_ext b, const uint idx)
{
    a->value[idx] = b;
}


coord_type get_coord_type_vector_item(const coord_type_vector a, const uint idx)
{
    coord_type tmp= a.value[idx];
    return tmp;
}


coord_type_ext get_coord_type_vector_ext_item(const coord_type_vector_ext a, const uint idx)
{
    coord_type_ext tmp= a.value[idx];
    return tmp;
}


void set_kd_tree_type_items(kdTree_type *a,
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


void get_kd_tree_type_items(kdTree_type a,
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

	*idx        = a.idx;
	*wgtCent 	= a.wgtCent;
	*midPoint 	= a.midPoint;
	*bnd_hi		= a.bnd_hi;
	*bnd_lo		= a.bnd_lo;
	*sum_sq 	= a.sum_sq;
	*count 		= a.count;
	*left 		= a.left;
	*right 		= a.right;

}

// conversion from data_type_ext to data_type
data_type conv_long_to_short(data_type_ext p)
{
    data_type result;
    for (uint d=0; d<D; d++) {
        coord_type tmp = (coord_type)get_coord_type_vector_ext_item(p.value,d);
        set_coord_type_vector_item(&result.value,tmp,d);
    }
    return result;
}

// conversion from data_type to data_type_ext
data_type_ext conv_short_to_long(data_type p)
{
    data_type_ext result;
    for (uint d=0; d<D; d++) {
        coord_type_ext tmp = (coord_type_ext)get_coord_type_vector_item(p.value,d);
        set_coord_type_vector_ext_item(&result.value,tmp,d);
    }
    return result;
}


void my_wait() {
	#ifdef INCLUDE_WAIT
	//ap_wait();
	#endif
}


void wait_for_n_cycles(const uint n)
{
	for (uint i=0; i<n; i++) {
		my_wait();
	}
}



mul_input_type saturate_mul_input(coord_type_ext val)
{
    if (val > MUL_MAX_VAL) {
        val = MUL_MAX_VAL;
    } else if (val < MUL_MIN_VAL) {
        val = MUL_MIN_VAL;
    }
    return (mul_input_type)val;
}


// fixed-point multiplication with saturation and scaling
coord_type_ext fi_mul(coord_type_ext op1, coord_type_ext op2)
{
    mul_input_type tmp_op1 = saturate_mul_input(op1);
    mul_input_type tmp_op2 = saturate_mul_input(op2);

    long int result_unscaled;
    result_unscaled = tmp_op1*tmp_op2;

    long int result_scaled;
    result_scaled = result_unscaled >> MUL_FRACTIONAL_BITS;
    coord_type_ext result;
    result = (coord_type_ext)result_scaled;
    return result;
}


// tree adder
coord_type_ext tree_adder(coord_type_ext *input_array)
{

    for(uint j=0;j<ceil(log2(D));j++)
    {
        if (j<ceil(log2(D))-1) {
            for(uint i = 0; i < uint(D/(1<<(j+1))); i++)
            {
                coord_type_ext tmp1 = input_array[2*i];
                coord_type_ext tmp2 = input_array[2*i+1];
                coord_type_ext tmp3 = tmp1+tmp2;
                input_array[i] = tmp3;
            }
            if (D > uint(D/(1<<(j+1)))*(1<<(j+1)) ) {
                input_array[uint(D/(1<<(j+1)))] = input_array[uint(D/(1<<(j+1))-1)*2+2];
            }
        }
        if (j== ceil(log2(D))-1) {
            coord_type_ext tmp1 = input_array[0];
            coord_type_ext tmp2 = input_array[1];
            coord_type_ext tmp3 = tmp1+tmp2;
            input_array[0] = tmp3;
        }
    }
    return input_array[0];
}

// tree adder (overloaded function)
coord_type_ext tree_adder(coord_type_ext *input_array,const uint m)
{

	for(uint j=0;j<ceil(log2(m));j++)
	{
			if (j<ceil(log2(m))-1) {
					for(uint i = 0; i < uint(m/(1<<(j+1))); i++)
					{
							coord_type_ext tmp1 = input_array[2*i];
							coord_type_ext tmp2 = input_array[2*i+1];
							coord_type_ext tmp3 = tmp1+tmp2;
							input_array[i] = tmp3;
					}
					if (m > uint(m/(1<<(j+1)))*(1<<(j+1)) ) {
							input_array[uint(m/(1<<(j+1)))] = input_array[uint(m/(1<<(j+1))-1)*2+2];
					}
			}
			if (j== ceil(log2(m))-1) {
					coord_type_ext tmp1 = input_array[0];
					coord_type_ext tmp2 = input_array[1];
					coord_type_ext tmp3 = tmp1+tmp2;
					input_array[0] = tmp3;
			}
	}
	return input_array[0];
}


// inner product of p1 and p2
void dot_product(data_type_ext p1,data_type_ext p2, coord_type_ext *r)
{

    coord_type_ext tmp_mul_res[D];

    for (uint d=0;d<D;d++) {
        mul_input_type tmp_op1 = saturate_mul_input(get_coord_type_vector_ext_item(p1.value,d));
        mul_input_type tmp_op2 = saturate_mul_input(get_coord_type_vector_ext_item(p2.value,d));
        coord_type_ext tmp_mul = tmp_op1*tmp_op2;
        tmp_mul_res[d] = tmp_mul;
    }

    *r = tree_adder(tmp_mul_res);
}


// compute the Euclidean distance between p1 and p2
void compute_distance(data_type_ext p1, data_type_ext p2, coord_type_ext *dist)
{

    data_type_ext tmp_p1 = p1;
    data_type_ext tmp_p2 = p2;
    coord_type_ext tmp_mul_res[D];

    for (uint d=0; d<D; d++) {
        coord_type_ext tmp_sub1 = get_coord_type_vector_ext_item(tmp_p1.value,d);
        coord_type_ext tmp_sub2 = get_coord_type_vector_ext_item(tmp_p2.value,d);
        coord_type_ext tmp = tmp_sub1 - tmp_sub2;
        coord_type_ext tmp_mul = fi_mul(tmp,tmp);
        tmp_mul_res[d] = tmp_mul;
    }

    *dist = tree_adder(tmp_mul_res);
}



// check whether any point of bounding box is closer to z than to z*
// this is a modified version of David Mount's code (http://www.cs.umd.edu/~mount/Projects/KMeans/ )
void tooFar_fi(data_type closest_cand, data_type cand, data_type bnd_lo, data_type bnd_hi, bool *too_far)
{

    coord_type_ext boxDot;
    coord_type_ext ccDot;

    data_type_ext tmp_closest_cand  = conv_short_to_long(closest_cand);
    data_type_ext tmp_cand          = conv_short_to_long(cand);
    data_type_ext tmp_bnd_lo        = conv_short_to_long(bnd_lo);
    data_type_ext tmp_bnd_hi        = conv_short_to_long(bnd_hi);

    coord_type_ext tmp_mul_res[D];
    coord_type_ext tmp_mul_res2[D];

    for (uint d = 0; d<D; d++) {
        coord_type_ext tmp_sub_op1 =  get_coord_type_vector_ext_item(tmp_cand.value,d);
        coord_type_ext tmp_sub_op2 =  get_coord_type_vector_ext_item(tmp_closest_cand.value,d);
        coord_type_ext ccComp = tmp_sub_op1-tmp_sub_op2;
        coord_type_ext tmp_mul = fi_mul(ccComp,ccComp);
        tmp_mul_res[d] = tmp_mul;

        coord_type_ext tmp_diff2;
        coord_type_ext tmp_sub_op3;
        coord_type_ext tmp_sub_op4;
        if (ccComp > 0) {
            tmp_sub_op3 = get_coord_type_vector_ext_item(tmp_bnd_hi.value,d);
        }
        else {
            tmp_sub_op3 = get_coord_type_vector_ext_item(tmp_bnd_lo.value,d);
        }
        tmp_sub_op4 = get_coord_type_vector_ext_item(tmp_closest_cand.value,d);
        tmp_diff2 = tmp_sub_op3 - tmp_sub_op4;

        coord_type_ext tmp_mul2 = fi_mul(tmp_diff2,ccComp);
        tmp_mul_res2[d] = tmp_mul2;

    }

    ccDot = tree_adder(tmp_mul_res);
    boxDot = tree_adder(tmp_mul_res2);

    coord_type_ext tmp_boxDot = boxDot<<1;
    bool tmp_res;
    if (ccDot>tmp_boxDot) {
        tmp_res = true;
    } else {
        tmp_res = false;
    }

    *too_far = tmp_res;
}

