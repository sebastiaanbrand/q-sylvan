#include <stdio.h>
#include <stdlib.h>
#include "sylvan_edge_weights_mpc.h"

#define UNUSED(x) (void)(x) // compile with unused params


/*****************<Implementation of edge_weights interface>*******************/

mpc_ptr
weight_mpc_malloc()
{
    // TODO
    return NULL;
}

void
_weight_mpc_value(void *wgt_store, EVBDD_WGT a, mpc_ptr res)
{
    UNUSED(wgt_store);
    UNUSED(a);
    UNUSED(res);
    // TODO
}

EVBDD_WGT
_weight_mpc_lookup_ptr(mpc_ptr a, void *wgt_store)
{
    UNUSED(a);
    UNUSED(wgt_store);
    // TODO
    return EVBDD_ZERO;
}

EVBDD_WGT
weight_mpc_lookup(mpc_ptr a)
{
    return _weight_mpc_lookup_ptr(a, wgt_storage);
}

void
init_mpc_one_zero(void *wgt_store)
{
    UNUSED(wgt_store);
    // TODO
}

void
weight_mpc_abs(mpc_ptr a)
{
    UNUSED(a);
    // TODO
}

void
weight_mpc_neg(mpc_ptr a)
{
    UNUSED(a);
    // TODO
}

void
weight_mpc_conj(mpc_ptr a)
{
    UNUSED(a);
    // TODO
}

void
weight_mpc_sqr(mpc_ptr a)
{
    weight_mpc_mul(a, a);
}

void
weight_mpc_add(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
}

void
weight_mpc_sub(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
}

void
weight_mpc_mul(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
}

void
weight_mpc_div(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
}

bool
weight_mpc_eq(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
    return false;
}

bool
weight_mpc_eps_close(mpc_ptr a, mpc_ptr b, double eps)
{
    UNUSED(a);
    UNUSED(b);
    UNUSED(eps);
    // TODO
    return false;
}

bool
weight_mpc_greater(mpc_ptr a, mpc_ptr b)
{
    UNUSED(a);
    UNUSED(b);
    // TODO
    return false;
}

EVBDD_WGT
wgt_mpc_norm_L2(EVBDD_WGT *low, EVBDD_WGT *high)
{
    UNUSED(low);
    UNUSED(high);
    // TODO
    return EVBDD_ZERO;
}

EVBDD_WGT
wgt_mpc_get_low_L2normed(EVBDD_WGT high)
{
    UNUSED(high);
    // TODO
    return EVBDD_ZERO;
}

void
weight_mpc_fprint(FILE *stream, mpc_ptr a)
{
    UNUSED(stream);
    UNUSED(a);
    // TODO
}

/*****************</Implementation of edge_weights interface>******************/
