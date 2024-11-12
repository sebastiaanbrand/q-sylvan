#ifndef WGT_MPC_H
#define WGT_MPC_H

#include "sylvan_edge_weights.h"
#include "edge_weight_storage/flt.h"


/******************<Implementation of edge_weights interface>******************/

mpc_ptr weight_mpc_malloc();
void _weight_mpc_value(void *wgt_store, EVBDD_WGT a, mpc_ptr res);
EVBDD_WGT weight_mpc_lookup(mpc_ptr a);
EVBDD_WGT _weight_mpc_lookup_ptr(mpc_ptr a, void *wgt_store);

void init_mpc_one_zero(void *wgt_store);

void weight_mpc_abs(mpc_ptr a);
void weight_mpc_neg(mpc_ptr a);
void weight_mpc_conj(mpc_ptr a);
void weight_mpc_sqr(mpc_ptr a);
void weight_mpc_add(mpc_ptr a, mpc_ptr b);
void weight_mpc_sub(mpc_ptr a, mpc_ptr b);
void weight_mpc_mul(mpc_ptr a, mpc_ptr b);
void weight_mpc_div(mpc_ptr a, mpc_ptr b);
bool weight_mpc_eq(mpc_ptr a, mpc_ptr b);
bool weight_mpc_eps_close(mpc_ptr a, mpc_ptr b, double eps);
bool weight_mpc_greater(mpc_ptr a, mpc_ptr b);

EVBDD_WGT wgt_mpc_norm_L2(EVBDD_WGT *low, EVBDD_WGT *high);
EVBDD_WGT wgt_mpc_get_low_L2normed(EVBDD_WGT high);

void weight_mpc_fprint(FILE *stream, mpc_ptr a);

/*****************</Implementation of edge_weights interface>******************/

#endif
