#include <stdio.h>

#include <sylvan_edge_weights.h>
#include <sylvan_edge_weights_complex.h>
#include <sylvan_int.h>


void *wgt_storage; // TODO: move to source file?
void *wgt_storage_new;

AADD_WGT AADD_ONE;
AADD_WGT AADD_ZERO;
AADD_WGT AADD_MIN_ONE;

/**********************<Managing the edge weight table>************************/

void sylvan_init_edge_weights(size_t size, double tol, edge_weight_type_t edge_weight_type, wgt_storage_backend_t backend);
void init_edge_weight_functions(edge_weight_type_t edge_weight_type);
void init_edge_weight_storage(size_t size, double tol, wgt_storage_backend_t backend, void **wgt_store);
void (*init_wgt_table_entries)(); // set by sylvan_init_aadd
uint64_t sylvan_get_edge_weight_table_size();
double sylvan_edge_weights_tolerance();
uint64_t sylvan_edge_weights_count_entries();
void sylvan_edge_weights_free();

/******************<Interface for different edge_weight_types>*****************/

weight_t (*weight_malloc)();
void (*_weight_value)(); // weight_value(WGT a, weight_type * res)
AADD_WGT (*weight_lookup)(); // weight_lookup_(weight_type a)
AADD_WGT (*_weight_lookup_ptr)(); // _weight_lookup_ptr(weight_type *a, void *wgt_store)
void (*init_one_zero)();

/* Arithmetic operations on edge weights */
void (*weight_abs)(); // a <-- |a|
void (*weight_neg)(); // a <-- -a
void (*weight_ccj)(); // a <-- complex conjugate of a
void (*weight_sqr)(); // a <-- a^2
void (*weight_add)(); // a <-- a + b
void (*weight_sub)(); // a <-- a - b
void (*weight_mul)(); // a <-- a * b
void (*weight_div)(); // a <-- a / b
bool (*weight_eq)(); // returns true iff a == b
bool (*weight_eps_close)(); // returns true iff dist(a,b) < eps
bool (*weight_greater)(); // returns true iff a > b

/* Normalization methods */
AADD_WGT (*wgt_norm_L2)(); // wgt_norm_L2(AADD_WGT *low, AADD_WGT *high)
AADD_WGT (*wgt_get_low_L2normed)(); // wgt_get_low_L2normed(AADD_WGT high)

void (*weight_fprint)(); // weight_fprint(FILE *stream, weight_t a)

/**********************<Managing the edge weight table>************************/

// Table parameters
static const double default_tolerance = 1e-14;
static double tolerance;
static wgt_storage_backend_t wgt_backend;
size_t table_size;

void sylvan_init_edge_weights(size_t size, double tol, edge_weight_type_t edge_weight_type, wgt_storage_backend_t backend)
{
    init_edge_weight_functions(edge_weight_type);
    init_edge_weight_storage(size, tol, backend, &wgt_storage);
    init_edge_weight_storage_gc();
}

void init_edge_weight_functions(edge_weight_type_t edge_weight_type)
{
    switch (edge_weight_type)
    {
    case WGT_COMPLEX_128:
        weight_malloc       = &weight_complex_malloc;
        _weight_value       = &_weight_complex_value;
        weight_lookup       = &weight_complex_lookup;
        _weight_lookup_ptr  = &_weight_complex_lookup_ptr;
        init_one_zero       = &init_complex_one_zero;
        weight_abs          = &weight_complex_abs;
        weight_neg          = &weight_complex_neg;
        weight_ccj          = &weight_complex_ccj;
        weight_sqr          = &weight_complex_sqr;
        weight_add          = &weight_complex_add;
        weight_sub          = &weight_complex_sub;
        weight_mul          = &weight_complex_mul;
        weight_div          = &weight_complex_div;
        weight_eq           = &weight_complex_eq;
        weight_eps_close    = &weight_complex_eps_close;
        weight_greater      = &weight_complex_greater;
        wgt_norm_L2         = &wgt_complex_norm_L2;
        wgt_get_low_L2normed= &wgt_complex_get_low_L2normed;
        weight_fprint       = &weight_complex_fprint;
        break;
    default:
        printf("ERROR: Unrecognized weight type = %d\n", edge_weight_type);
        exit(1);
        break;
    }
    
}

void
init_edge_weight_storage(size_t size, double tol, wgt_storage_backend_t backend, void **wgt_store)
{
    tolerance = (tol < 0) ? default_tolerance : tol;
    table_size = size;
    wgt_backend = backend;

    init_wgt_storage_functions(backend);

    // create actual table
    *wgt_store = wgt_store_create(table_size, tolerance);

    // Set AADD_WGT values for 1, 0 (and -1)
    init_one_zero(*wgt_store);
}

uint64_t
sylvan_get_edge_weight_table_size()
{
    return table_size;
}

double
sylvan_edge_weights_tolerance()
{
    return wgt_store_get_tol();
}

uint64_t
sylvan_edge_weights_count_entries()
{
    return wgt_store_num_entries(wgt_storage);
}

void
sylvan_edge_weights_free()
{
    wgt_store_free(wgt_storage);
}

/*********************</Managing the edge weight table>************************/





/*************************<GC of edge weight table>****************************/

// Keep estimate for number of entries for gc purposes
size_t table_entries_est;
DECLARE_THREAD_LOCAL(table_entries_local, size_t); // these are added to _est
static const uint64_t table_entries_local_buffer = 1000; // every 1000 entries

void
init_edge_weight_storage_gc()
{
    LOCALIZE_THREAD_LOCAL(table_entries_local, size_t);
    // NOTE: the sum of the local counters sometimes exceeds the actual total
    // number of entries (when just counting the global value with atomic adds
    // after every insert). This might be because 'ctable_entries_local = 0' 
    // doesn't set it to 0 for all threads. Since it's only an estimate it is
    // not a huge issue though.
    // TODO: figure out how to get lace to handle this counting better
    table_entries_est = 0;
    table_entries_local = 0;
    (void) table_entries_local;
}

uint64_t
wgt_table_entries_estimate()
{
    return table_entries_est;
}

void wgt_table_gc_inc_entries_estimate()
{
    LOCALIZE_THREAD_LOCAL(table_entries_local, size_t);
    table_entries_local += 1;
    if (table_entries_local >= table_entries_local_buffer) {
        __sync_fetch_and_add(&table_entries_est, table_entries_local);
        table_entries_local = 0;
    }
}

void
wgt_table_gc_init_new(void (*init_wgt_table_entries)())
{
    // init new (empty) edge weight storage
    init_edge_weight_storage(table_size, tolerance, wgt_backend, &wgt_storage_new);

    // reset estimate entries counters
    LOCALIZE_THREAD_LOCAL(table_entries_local, size_t);
    table_entries_est = 0;
    table_entries_local = 0;
    (void) table_entries_local;

    // Fill new with initial values (temp rename to wgt_store because
    // init_wgt_table_entries initializes wgt_storage, not wgt_storage_new)
    void *wgt_store_tmp = wgt_storage;
    wgt_storage = wgt_storage_new;
    if (init_wgt_table_entries != NULL) {
        init_wgt_table_entries();
    }
    wgt_storage_new = wgt_storage;  // wgt_store_new now has initial values
    wgt_storage = wgt_store_tmp;    // wgt_storage now has the old values
}

void
wgt_table_gc_delete_old()
{
    // delete  old (full) table + set new as current
    wgt_store_free(wgt_storage);
    wgt_storage = wgt_storage_new;
}

AADD_WGT
wgt_table_gc_keep(AADD_WGT a)
{
    // move from current (old) to new
    weight_t wa = weight_malloc();
    _weight_value(wgt_storage, a, wa);
    AADD_WGT res = _weight_lookup_ptr(wa, wgt_storage_new);
    free(wa);
    return res;
}

/************************</GC of edge weight table>****************************/





/********************<For caching arithmetic operations>***********************/

static bool CACHE_WGT_OPS = true;
static bool CACHE_INV_OPS = true;

void
wgt_set_inverse_chaching(bool on)
{
    CACHE_INV_OPS = on;
}

static void
order_inputs(AADD_WGT *a, AADD_WGT *b) 
{
    AADD_WGT x = (*a > *b) ? *a : *b;
    AADD_WGT y = (*a > *b) ? *b : *a;
    *a = x;
    *b = y;
}

static bool
cache_get_add(AADD_WGT a, AADD_WGT b, AADD_WGT *res)
{
    order_inputs(&a, &b);
    if (cache_get3(CACHE_WGT_ADD, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_ADD_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_add(AADD_WGT a, AADD_WGT b, AADD_WGT res)
{
    order_inputs(&a, &b);
    if (cache_put3(CACHE_WGT_ADD, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_ADD_CACHEDPUT);
    }
}

static void
cache_put_sub(AADD_WGT a, AADD_WGT b, AADD_WGT res)
{
    if (cache_put3(CACHE_WGT_SUB, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_SUB_CACHEDPUT);
    }
}

static bool
cache_get_sub(AADD_WGT a, AADD_WGT b, AADD_WGT *res)
{
    if (cache_get3(CACHE_WGT_SUB, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_SUB_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_mul(AADD_WGT a, AADD_WGT b, AADD_WGT res)
{
    order_inputs(&a, &b);
    if (cache_put3(CACHE_WGT_MUL, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_MUL_CACHEDPUT);
    }
    if (CACHE_INV_OPS) {
        // put inverse as well (empirically seems not so beneficial)
        if (cache_put3(CACHE_WGT_DIV, res, b, sylvan_false, a)) {
            sylvan_stats_count(WGT_DIV_CACHEDPUT);
        }
        if (cache_put3(CACHE_WGT_DIV, res, a, sylvan_false, b)) {
            sylvan_stats_count(WGT_DIV_CACHEDPUT);
        }
    }
}

static bool
cache_get_mul(AADD_WGT a, AADD_WGT b, AADD_WGT *res)
{
    order_inputs(&a, &b);
    if (cache_get3(CACHE_WGT_MUL, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_MUL_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_div(AADD_WGT a, AADD_WGT b, AADD_WGT res)
{
    if (cache_put3(CACHE_WGT_DIV, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_DIV_CACHEDPUT);
    }
    if (CACHE_INV_OPS) {
        // put inverse as well (empirically seems beneficial)
        order_inputs(&b, &res);
        if (cache_put3(CACHE_WGT_MUL, b, res, sylvan_false, a)) {
            sylvan_stats_count(WGT_MUL_CACHEDPUT);
        }
    }
}

static bool
cache_get_div(AADD_WGT a, AADD_WGT b, AADD_WGT *res)
{
    if (cache_get3(CACHE_WGT_DIV, a, b, sylvan_false, res)) {
        sylvan_stats_count(WGT_DIV_CACHED);
        return true;
    }
    return false;
}

/*******************</For caching arithmetic operations>***********************/





/*********************<Arithmetic functions on AADD_WGT's>*********************/

AADD_WGT
wgt_abs(AADD_WGT a)
{
    // special cases
    if (a == AADD_ZERO || a == AADD_ONE) return a;
    if (a == AADD_MIN_ONE) return AADD_ONE;

    AADD_WGT res;

    weight_t w = weight_malloc();
    weight_value(a, w);
    weight_abs(w);
    res = weight_lookup_ptr(w);
    free(w);

    return res;
}

AADD_WGT 
wgt_neg(AADD_WGT a)
{
    // special cases
    if (a == AADD_ZERO) return AADD_ZERO;
    if (a == AADD_ONE) return AADD_MIN_ONE;
    if (a == AADD_MIN_ONE) return AADD_ONE;

    AADD_WGT res;

    weight_t w = weight_malloc();
    weight_value(a, w);
    weight_neg(w);
    res = weight_lookup_ptr(w);
    free(w);

    return res; 
}

AADD_WGT
wgt_ccj(AADD_WGT a)
{
    // special cases
    if (a == AADD_ZERO || a == AADD_ONE || a == AADD_MIN_ONE) return a;

    AADD_WGT res;

    weight_t w = weight_malloc();
    weight_value(a, w);
    weight_ccj(w);
    res = weight_lookup_ptr(w);
    free(w);

    return res;
}

AADD_WGT
wgt_add(AADD_WGT a, AADD_WGT b)
{
    // special cases
    if (a == AADD_ZERO) return b;
    if (b == AADD_ZERO) return a;

    // check cache
    AADD_WGT res;
    if (CACHE_WGT_OPS) {
        if (cache_get_add(a, b, &res)) return res;
    }

    // compute and lookup result in edge weight table
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();
    weight_value(a, wa);
    weight_value(b, wb);
    weight_add(wa, wb);
    res = weight_lookup_ptr(wa);
    free(wa);
    free(wb);

    // insert in cache
    if (CACHE_WGT_OPS) {
        cache_put_add(a, b, res);
    }
    return res;
}

AADD_WGT
wgt_sub(AADD_WGT a, AADD_WGT b)
{
    // special cases
    if (b == AADD_ZERO) return a;
    if (a == AADD_ZERO) return wgt_neg(b);

    // check cache
    AADD_WGT res;
    if (CACHE_WGT_OPS) {
        if (cache_get_sub(a, b, &res)) return res;
    }

    // compute and lookup result in edge weight table
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();
    weight_value(a, wa);
    weight_value(b, wb);
    weight_sub(wa, wb);
    res = weight_lookup_ptr(wa);
    free(wa);
    free(wb);

    // insert in cache
    if (CACHE_WGT_OPS) {
        cache_put_sub(a, b, res);
    }
    return res;
}

AADD_WGT
wgt_mul(AADD_WGT a, AADD_WGT b)
{
    // special cases
    if (a == AADD_ONE) return b;
    if (b == AADD_ONE) return a;
    if (a == AADD_ZERO || b == AADD_ZERO) return AADD_ZERO;

    // check cache
    AADD_WGT res;
    if (CACHE_WGT_OPS) {
        if (cache_get_mul(a, b, &res)) return res;
    }

    // compute and lookup result in edge weight table
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();
    weight_value(a, wa);
    weight_value(b, wb);
    weight_mul(wa, wb);
    res = weight_lookup_ptr(wa);
    free(wa);
    free(wb);

    // insert in cache
    if (CACHE_WGT_OPS) {
        cache_put_mul(a, b, res);
    }
    return res;
}

AADD_WGT
wgt_div(AADD_WGT a, AADD_WGT b)
{
    // special cases
    if (a == b)         return AADD_ONE;
    if (a == AADD_ZERO) return AADD_ZERO;
    if (b == AADD_ONE)  return a;

    // check cache
    AADD_WGT res;
    if (CACHE_WGT_OPS) {
        if (cache_get_div(a, b, &res)) return res;
    }

    // compute and lookup result in edge weight table
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();
    weight_value(a, wa);
    weight_value(b, wb);
    weight_div(wa, wb);
    res = weight_lookup_ptr(wa);
    free(wa);
    free(wb);

    // insert in cache
    if (CACHE_WGT_OPS) {
        cache_put_div(a, b, res);
    }
    return res;
}

/********************</Arithmetic functions on AADD_WGT's>*********************/





/*************************<Comparators on AADD_WGT's>**************************/

bool
wgt_eq(AADD_WGT a, AADD_WGT b)
{
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();

    weight_value(a, wa);
    weight_value(b, wb);
    bool res = weight_eq(wa, wb);

    free(wa);
    free(wb);

    return res;
}

bool
wgt_eps_close(AADD_WGT a, AADD_WGT b, double eps)
{
    weight_t wa = weight_malloc();
    weight_t wb = weight_malloc();

    weight_value(a, wa);
    weight_value(b, wb);
    bool res = weight_eps_close(wa, wb, eps);

    free(wa);
    free(wb);

    return res;
}

bool
wgt_approx_eq(AADD_WGT a, AADD_WGT b)
{
    return wgt_eps_close(a, b, wgt_store_get_tol());
}

/************************</Comparators on AADD_WGT's>**************************/





/*************************<Edge weight normalization>**************************/

AADD_WGT
wgt_norm_low(AADD_WGT *low, AADD_WGT *high)
{
    // Normalize using low if low != 0
    AADD_WGT norm;
    if(*low != AADD_ZERO){
        *high = wgt_div(*high, *low);
        norm  = *low;
        *low  = AADD_ONE;
    }
    else {
        norm  = *high;
        *high = AADD_ONE;
    }
    return norm;    
}

AADD_WGT
wgt_norm_largest(AADD_WGT *low, AADD_WGT *high)
{
    AADD_WGT norm;
    if (*low == *high) {
        norm  = *low;
        *low  = AADD_ONE;
        *high = AADD_ONE;
        return norm;
    }

    // Normalize using the absolute greatest value
    weight_t wl = weight_malloc();
    weight_t wh = weight_malloc();
    weight_value(*low,  wl);
    weight_value(*high, wh);

    if (weight_greater(wh, wl)) {
        // high greater than low, divide both by high
        *low = wgt_div(*low, *high);
        norm  = *high;
        *high = AADD_ONE;
    }
    else {
        // low greater than high (or equal magnitude), divide both by low
        *high = wgt_div(*high, *low);
        norm = *low;
        *low  = AADD_ONE;
    }

    free(wl);
    free(wh);
    return norm;
}

/************************</Edge weight normalization>**************************/





/************************<Printing & utility functions>************************/

void wgt_fprint(FILE *stream, AADD_WGT a)
{
    weight_t w = weight_malloc();
    weight_value(a, w);
    weight_fprint(stream, w);
    free(w);
}

/************************<Printing & utility functions>************************/
