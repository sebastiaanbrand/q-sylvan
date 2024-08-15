/*
 * Copyright 2011-2016 Formal Methods and Tools, University of Twente
 * Copyright 2016-2017 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <sylvan_int.h>
#include <sylvan_aadd.h>
#include <sylvan_refs.h>

static int granularity = 1; // operation cache access granularity


bool larger_wgt_indices;
int weight_norm_strat;
AADD_WGT (*normalize_weights)(AADD_WGT *, AADD_WGT *);

/*******************<Garbage collection, references, marking>******************/

/**
 * Most of this gc code is copy-paste from sylvan_mtbdd.c, however because the 
 * bit structure of BDDs and AADDs are a bit different we can't use the mtbdd
 * code directly. Since the way sylvan_mtbdd/sylvan_common/sylvan_table are
 * structured we need to copy-paste a few more functions/variables than we 
 * actually change.
 */

/* 
 * Recursively mark AADD nodes as 'in use'.
 * This is really the only gc function which is different for AADDs vs MTBDDs.
 */
VOID_TASK_IMPL_1(aadd_gc_mark_rec, AADD, a)
{
    if (AADD_TARGET(a) == AADD_TERMINAL) return;

    if (llmsset_mark(nodes, AADD_TARGET(a))) {
        aaddnode_t n = AADD_GETNODE(AADD_TARGET(a));
        SPAWN(aadd_gc_mark_rec, aaddnode_getptrlow(n));
        CALL(aadd_gc_mark_rec, aaddnode_getptrhigh(n));
        SYNC(aadd_gc_mark_rec);
    }
}

/**
 * External references
 */
refs_table_t aadd_refs;
refs_table_t aadd_protected;
static int aadd_protected_created = 0;

void
aadd_protect(AADD *a)
{
    if (!aadd_protected_created) {
        // In C++, sometimes mtbdd_protect is called before Sylvan is initialized. Just create a table.
        protect_create(&aadd_protected, 4096);
        aadd_protected_created = 1;
    }
    protect_up(&aadd_protected, (size_t)a);
}

void
aadd_unprotect(AADD *a)
{
    if (aadd_protected.refs_table != NULL) protect_down(&aadd_protected, (size_t)a);
}

size_t
aadd_count_protected()
{
    return protect_count(&aadd_protected);
}

/* Called during garbage collection */
VOID_TASK_0(aadd_gc_mark_external_refs)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = refs_iter(&aadd_refs, 0, aadd_refs.refs_size);
    while (it != NULL) {
        SPAWN(aadd_gc_mark_rec, refs_next(&aadd_refs, &it, aadd_refs.refs_size));
        count++;
    }
    while (count--) {
        SYNC(aadd_gc_mark_rec);
    }
}

/* Called during garbage collection */
VOID_TASK_0(aadd_gc_mark_protected)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = protect_iter(&aadd_protected, 0, aadd_protected.refs_size);
    while (it != NULL) {
        AADD *to_mark = (AADD*)protect_next(&aadd_protected, &it, aadd_protected.refs_size);
        SPAWN(aadd_gc_mark_rec, *to_mark);
        count++;
    }
    while (count--) {
        SYNC(aadd_gc_mark_rec);
    }
}

/* Infrastructure for internal markings */
typedef struct aadd_refs_task
{
    Task *t;
    void *f;
} *aadd_refs_task_t;

typedef struct aadd_refs_internal
{
    const AADD **pbegin, **pend, **pcur;
    AADD *rbegin, *rend, *rcur;
    aadd_refs_task_t sbegin, send, scur;
} *aadd_refs_internal_t;

DECLARE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);

VOID_TASK_2(aadd_refs_mark_p_par, const AADD**, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            aadd_gc_mark_rec(**(begin++));
            count--;
        }
    } else {
        SPAWN(aadd_refs_mark_p_par, begin, count / 2);
        CALL(aadd_refs_mark_p_par, begin + (count / 2), count - count / 2);
        SYNC(aadd_refs_mark_p_par);
    }
}

VOID_TASK_2(aadd_refs_mark_r_par, AADD*, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            aadd_gc_mark_rec(*begin++);
            count--;
        }
    } else {
        SPAWN(aadd_refs_mark_r_par, begin, count / 2);
        CALL(aadd_refs_mark_r_par, begin + (count / 2), count - count / 2);
        SYNC(aadd_refs_mark_r_par);
    }
}

VOID_TASK_2(aadd_refs_mark_s_par, aadd_refs_task_t, begin, size_t, count)
{
    if (count < 32) {
        while (count > 0) {
            Task *t = begin->t;
            if (!TASK_IS_STOLEN(t)) return;
            if (t->f == begin->f && TASK_IS_COMPLETED(t)) {
                aadd_gc_mark_rec(*(AADD*)TASK_RESULT(t));
            }
            begin += 1;
            count -= 1;
        }
    } else {
        if (!TASK_IS_STOLEN(begin->t)) return;
        SPAWN(aadd_refs_mark_s_par, begin, count / 2);
        CALL(aadd_refs_mark_s_par, begin + (count / 2), count - count / 2);
        SYNC(aadd_refs_mark_s_par);
    }
}

VOID_TASK_0(aadd_refs_mark_task)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    SPAWN(aadd_refs_mark_p_par, aadd_refs_key->pbegin, aadd_refs_key->pcur-aadd_refs_key->pbegin);
    SPAWN(aadd_refs_mark_r_par, aadd_refs_key->rbegin, aadd_refs_key->rcur-aadd_refs_key->rbegin);
    CALL(aadd_refs_mark_s_par, aadd_refs_key->sbegin, aadd_refs_key->scur-aadd_refs_key->sbegin);
    SYNC(aadd_refs_mark_r_par);
    SYNC(aadd_refs_mark_p_par);
}

/* Called during garbage collection */
VOID_TASK_0(aadd_refs_mark)
{
    TOGETHER(aadd_refs_mark_task);
}

VOID_TASK_0(aadd_refs_init_task)
{
    aadd_refs_internal_t s = (aadd_refs_internal_t)malloc(sizeof(struct aadd_refs_internal));
    s->pcur = s->pbegin = (const AADD**)malloc(sizeof(AADD*) * 1024);
    s->pend = s->pbegin + 1024;
    s->rcur = s->rbegin = (AADD*)malloc(sizeof(AADD) * 1024);
    s->rend = s->rbegin + 1024;
    s->scur = s->sbegin = (aadd_refs_task_t)malloc(sizeof(struct aadd_refs_task) * 1024);
    s->send = s->sbegin + 1024;
    SET_THREAD_LOCAL(aadd_refs_key, s);
}

VOID_TASK_0(aadd_refs_init)
{
    INIT_THREAD_LOCAL(aadd_refs_key);
    TOGETHER(aadd_refs_init_task);
    sylvan_gc_add_mark(TASK(aadd_refs_mark));
}

VOID_TASK_0(aadd_refs_cleanup_task)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    free(aadd_refs_key->pbegin);
    free(aadd_refs_key->rbegin);
    free(aadd_refs_key->sbegin);
    free(aadd_refs_key);
}

/**
 * Called by aadd_quit. Cleans up the (thread local) malloc'ed aadd_refs_key.
 * 
 * NOTE: this cleanup isn't done in sylvan_mtbdd.c, but not doing this causes
 * memory leaks when calling initializing and quiting Sylvan multiple times
 * during the same program run.
 */
VOID_TASK_0(aadd_refs_cleanup)
{
    TOGETHER(aadd_refs_cleanup_task);
}

void
aadd_refs_ptrs_up(aadd_refs_internal_t aadd_refs_key)
{
    size_t cur = aadd_refs_key->pcur - aadd_refs_key->pbegin;
    size_t size = aadd_refs_key->pend - aadd_refs_key->pbegin;
    aadd_refs_key->pbegin = (const AADD**)realloc(aadd_refs_key->pbegin, sizeof(AADD*) * size * 2);
    aadd_refs_key->pcur = aadd_refs_key->pbegin + cur;
    aadd_refs_key->pend = aadd_refs_key->pbegin + (size * 2);
}

AADD __attribute__((noinline))
aadd_refs_refs_up(aadd_refs_internal_t aadd_refs_key, AADD res)
{
    long size = aadd_refs_key->rend - aadd_refs_key->rbegin;
    aadd_refs_key->rbegin = (AADD*)realloc(aadd_refs_key->rbegin, sizeof(AADD) * size * 2);
    aadd_refs_key->rcur = aadd_refs_key->rbegin + size;
    aadd_refs_key->rend = aadd_refs_key->rbegin + (size * 2);
    return res;
}

void __attribute__((noinline))
aadd_refs_tasks_up(aadd_refs_internal_t aadd_refs_key)
{
    long size = aadd_refs_key->send - aadd_refs_key->sbegin;
    aadd_refs_key->sbegin = (aadd_refs_task_t)realloc(aadd_refs_key->sbegin, sizeof(struct aadd_refs_task) * size * 2);
    aadd_refs_key->scur = aadd_refs_key->sbegin + size;
    aadd_refs_key->send = aadd_refs_key->sbegin + (size * 2);
}

void __attribute__((unused))
aadd_refs_pushptr(const AADD *ptr)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    *aadd_refs_key->pcur++ = ptr;
    if (aadd_refs_key->pcur == aadd_refs_key->pend) aadd_refs_ptrs_up(aadd_refs_key);
}

void __attribute__((unused))
aadd_refs_popptr(size_t amount)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    aadd_refs_key->pcur -= amount;
}

AADD __attribute__((unused))
aadd_refs_push(AADD a)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    *(aadd_refs_key->rcur++) = a;
    if (aadd_refs_key->rcur == aadd_refs_key->rend) return aadd_refs_refs_up(aadd_refs_key, a);
    else return a;
}

void __attribute__((unused))
aadd_refs_pop(long amount)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    aadd_refs_key->rcur -= amount;
}

void
aadd_refs_spawn(Task *t)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    aadd_refs_key->scur->t = t;
    aadd_refs_key->scur->f = t->f;
    aadd_refs_key->scur += 1;
    if (aadd_refs_key->scur == aadd_refs_key->send) aadd_refs_tasks_up(aadd_refs_key);
}

AADD
aadd_refs_sync(AADD result)
{
    LOCALIZE_THREAD_LOCAL(aadd_refs_key, aadd_refs_internal_t);
    aadd_refs_key->scur -= 1;
    return result;
}

/******************</Garbage collection, references, marking>******************/





/*************************<Cleaning edge weight table>*************************/

static int auto_gc_wgt_table  = 1;
static double wgt_table_gc_thres = 0.5;

void
aadd_set_auto_gc_wgt_table(bool enabled)
{
    auto_gc_wgt_table = enabled;
}

void
aadd_set_gc_wgt_table_thres(double fraction_filled)
{
    wgt_table_gc_thres = fraction_filled;
}

double
aadd_get_gc_wgt_table_thres()
{
    return wgt_table_gc_thres;
}


void
aadd_gc_wgt_table()
{
    // gc edge weight table and keep wgts of protected AADDs (and update those)
    // 1. Create new edge weight table table
    wgt_table_gc_init_new(init_wgt_table_entries);

    // 2. Fill new table with wgts in protected AADDs and update those AADDs
    uint64_t *it = protect_iter(&aadd_protected, 0, aadd_protected.refs_size);
    while (it != NULL) {
        AADD *to_protect_wgts = (AADD*)protect_next(&aadd_protected, &it, aadd_protected.refs_size);
        if (to_protect_wgts != NULL) {
            *to_protect_wgts = _fill_new_wgt_table(*to_protect_wgts);
        }
    }

    // 3. Delete old edge weight table
    wgt_table_gc_delete_old();

    // 4. Any cache we migh have is now invalid because the same edge weights 
    //    might now have different indices in the edge weight table
    sylvan_clear_cache();
}

TASK_IMPL_1(AADD, _fill_new_wgt_table, AADD, a)
{
    // Check cache
    AADD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_AADD_CLEAN_WGT_TABLE, 0LL, a, 0LL, &res)) {
            return res;
        }
    }

    // Move weight from old to new table, get new index
    AADD_WGT new_wgt = wgt_table_gc_keep(AADD_WEIGHT(a));
    a = aadd_bundle(AADD_TARGET(a), new_wgt);

    // If terminal, return
    if (AADD_TARGET(a) == AADD_TERMINAL) return a;
    
    // Recursive for children
    AADD low, high;
    aaddnode_t n = AADD_GETNODE(AADD_TARGET(a));
    aaddnode_getchilderen(n, &low, &high);
    aadd_refs_spawn(SPAWN(_fill_new_wgt_table, high));
    low = CALL(_fill_new_wgt_table, low);
    aadd_refs_push(low);
    high = aadd_refs_sync(SYNC(_fill_new_wgt_table));
    aadd_refs_pop(1);

    // We don't need to use the 'aadd_makenode()' function which normalizes the 
    // weights, because the AADD doesn't actually change, only the WGT indices,
    // but none of the actual values.
    AADD_TARG ptr = _aadd_makenode(aaddnode_getvar(n), AADD_TARGET(low), AADD_TARGET(high), AADD_WEIGHT(low), AADD_WEIGHT(high));

    // Put in cache, return
    res = aadd_bundle(ptr, new_wgt);
    if (cachenow) cache_put3(CACHE_AADD_CLEAN_WGT_TABLE, 0LL, a, 0LL, res);
    return res;
}

bool
aadd_test_gc_wgt_table()
{
    uint64_t entries = wgt_table_entries_estimate();
    uint64_t size    = sylvan_get_edge_weight_table_size();
    return ( ((double)entries / (double)size) > wgt_table_gc_thres );
}

/************************</Cleaning edge weight table>*************************/





/******************************<Initialization>********************************/

/**
 * Initialize and quit functions
 */
static int aadd_initialized = 0;

static void
aadd_quit()
{
    refs_free(&aadd_refs);
    if (aadd_protected_created) {
        protect_free(&aadd_protected);
        aadd_protected_created = 0;
    }
    RUN(aadd_refs_cleanup);
    aadd_initialized = 0;
    sylvan_edge_weights_free();
}

void
sylvan_init_aadd(size_t min_wgt_tablesize, size_t max_wgt_tablesize,
                 double wgt_tab_tolerance, int edge_weigth_backend, 
                 int norm_strat, void *init_wgt_tab_entries)
{
    if (aadd_initialized) return;
    aadd_initialized = 1;

    int index_size = (int) ceil(log2(max_wgt_tablesize));
    if (edge_weigth_backend == REAL_TUPLES_HASHMAP || edge_weigth_backend == REAL_TREE)
        index_size = index_size*2;
    if (index_size > 33) {
        fprintf(stderr,"max edge weight storage size is 2^33 (2^16 when using storing r and i seperately)\n");
        exit(1);
    }
    if (index_size > 23) larger_wgt_indices = true;
    else larger_wgt_indices = false;

    sylvan_register_quit(aadd_quit);
    sylvan_gc_add_mark(TASK(aadd_gc_mark_external_refs));
    sylvan_gc_add_mark(TASK(aadd_gc_mark_protected));

    refs_create(&aadd_refs, 1024);
    if (!aadd_protected_created) {
        protect_create(&aadd_protected, 4096);
        aadd_protected_created = 1;
    }

    // TODO: pass edge weight type to sylvan_init_aadd
    if (min_wgt_tablesize > max_wgt_tablesize) min_wgt_tablesize = max_wgt_tablesize;
    sylvan_init_edge_weights(min_wgt_tablesize, max_wgt_tablesize, 
                             wgt_tab_tolerance, WGT_COMPLEX_128, 
                             edge_weigth_backend);
    
    init_wgt_table_entries = init_wgt_tab_entries;
    if (init_wgt_table_entries != NULL) {
        init_wgt_table_entries();
    }

    weight_norm_strat = norm_strat;
    switch (norm_strat)
    {
    case NORM_LOW:
        normalize_weights = &wgt_norm_low;
        break;
    case NORM_MAX:
        normalize_weights = &wgt_norm_max;
        break;
    case NORM_MIN:
        normalize_weights = &wgt_norm_min;
        break;
    case NORM_L2:
        normalize_weights = wgt_norm_L2;
        break;
    default:
        printf("Edge weight normalization strategy not recognized\n");
        exit(1);
        break;
    }

    RUN(aadd_refs_init);
}

void
sylvan_init_aadd_defaults(size_t min_wgt_tablesize, size_t max_wgt_tablesize)
{
    sylvan_init_aadd(min_wgt_tablesize, max_wgt_tablesize, -1, COMP_HASHMAP, NORM_LOW, NULL);
}


void
aadd_set_caching_granularity(int g)
{
    granularity = g;
}

/*****************************</Initialization>********************************/





/**************************<Matrix/vector operations>**************************/

AADD_WGT
aadd_getvalue(AADD a, bool* path)
{
    AADD_WGT res = AADD_ONE;
    AADD low, high;
    for (;;) {
        res = wgt_mul(res, AADD_WEIGHT(a));
        
        // if the current edge is pointing to the terminal, we're done.
        if (AADD_TARGET(a) == AADD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
        BDDVAR var     = aaddnode_getvar(node);
        aaddnode_getchilderen(node, &low, &high);

        // Condition low/high choice on basis state vector[var]
        a = (path[var] == 0) ? low : high;
    }

    return res;
}

static void
aadd_do_before_mult(AADD *a, AADD *b)
{
    // check if edge weight table needs gc
    if (auto_gc_wgt_table && aadd_test_gc_wgt_table()) {
        aadd_protect(a);
        aadd_protect(b);
        aadd_gc_wgt_table();
        aadd_unprotect(a);
        aadd_unprotect(b);
    }
}

static void
norm_commuting_cache_key(AADD a, AADD b, AADD *x, AADD *y)
{
    if (a < b) {
        *x = a;
        *y = b;
    }
    else {
        *x = b;
        *y = a;
    }
}

TASK_IMPL_2(AADD, aadd_plus, AADD, a, AADD, b)
{
    // Trivial cases
    if(AADD_WEIGHT(a) == AADD_ZERO) return b;
    if(AADD_WEIGHT(b) == AADD_ZERO) return a;

    sylvan_gc_test();

    // Get var(a) and var(b)
    AADD low_a, low_b, high_a, high_b, res;
    BDDVAR var_a = UINT32_MAX, var_b = UINT32_MAX, topvar;
    if (AADD_TARGET(a) != AADD_TERMINAL) {
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
        var_a  = aaddnode_getvar(node);
    }
    if (AADD_TARGET(b) != AADD_TERMINAL) {
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(b));
        var_b  = aaddnode_getvar(node);
    }

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    aadd_get_topvar(a, var_b, &topvar, &low_a, &high_a);
    aadd_get_topvar(b, var_a, &topvar, &low_b, &high_b);

    // Base/terminal case: same target and same variable
    if(AADD_TARGET(a) == AADD_TARGET(b) && var_a == var_b){
        AADD_WGT sum = wgt_add(AADD_WEIGHT(a), AADD_WEIGHT(b));
        res = aadd_bundle(AADD_TARGET(a), sum);
        return res;
    }

    // Check cache
    AADD x, y;
    norm_commuting_cache_key(a, b, &x, &y); // (a+b) = (b+a) so normalize cache key
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_AADD_PLUS, sylvan_false, x, y, &res)) {
            sylvan_stats_count(AADD_PLUS_CACHED);
            return res;
        }
    }

    // If not base/terminal case, pass edge weight of current edge down
    AADD_WGT wgt_la, wgt_ha, wgt_lb, wgt_hb;
    wgt_la = wgt_mul(AADD_WEIGHT(a), AADD_WEIGHT(low_a));
    wgt_ha = wgt_mul(AADD_WEIGHT(a), AADD_WEIGHT(high_a));
    wgt_lb = wgt_mul(AADD_WEIGHT(b), AADD_WEIGHT(low_b));
    wgt_hb = wgt_mul(AADD_WEIGHT(b), AADD_WEIGHT(high_b));
    low_a  = aadd_refs_push(aadd_bundle(AADD_TARGET(low_a),  wgt_la));
    high_a = aadd_refs_push(aadd_bundle(AADD_TARGET(high_a), wgt_ha));
    low_b  = aadd_refs_push(aadd_bundle(AADD_TARGET(low_b),  wgt_lb));
    high_b = aadd_refs_push(aadd_bundle(AADD_TARGET(high_b), wgt_hb));

    // Recursive calls down
    aadd_refs_spawn(SPAWN(aadd_plus, high_a, high_b));
    AADD low = aadd_refs_push(CALL(aadd_plus, low_a, low_b));
    AADD high = aadd_refs_sync(SYNC(aadd_plus));
    aadd_refs_pop(5);

    // Put in cache, return
    res = aadd_makenode(topvar, low, high);
    if (cachenow) {
        if (cache_put3(CACHE_AADD_PLUS, sylvan_false, x, y, res)) 
            sylvan_stats_count(AADD_PLUS_CACHEDPUT);
    }
    return res;
}



/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(AADD, aadd_matvec_mult, AADD, mat, AADD, vec, BDDVAR, nvars)
{
    aadd_do_before_mult(&mat, &vec);
    aadd_refs_push(mat); aadd_refs_push(vec);
    AADD res = CALL(aadd_matvec_mult_rec, mat, vec, nvars, 0);
    aadd_refs_pop(2);
    return res;
}

/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(AADD, aadd_matmat_mult, AADD, a, AADD, b, BDDVAR, nvars)
{
    aadd_do_before_mult(&a, &b);
    aadd_refs_push(a); aadd_refs_push(b);
    AADD res = CALL(aadd_matmat_mult_rec, a, b, nvars, 0);
    aadd_refs_pop(2);
    return res;
}

TASK_IMPL_4(AADD, aadd_matvec_mult_rec, AADD, mat, AADD, vec, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (AADD_WEIGHT(mat) == AADD_ZERO || AADD_WEIGHT(vec) == AADD_ZERO)
        return aadd_bundle(AADD_TERMINAL, AADD_ZERO);
    
    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(AADD_TARGET(mat) == AADD_TERMINAL);
        assert(AADD_TARGET(vec) == AADD_TERMINAL);
        AADD_WGT prod = wgt_mul(AADD_WEIGHT(mat), AADD_WEIGHT(vec));
        return aadd_bundle(AADD_TERMINAL, prod);
    }

    sylvan_gc_test();

    // Check cache
    AADD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_AADD_MATVEC_MULT, nextvar, AADD_TARGET(mat), AADD_TARGET(vec), &res)) {
            sylvan_stats_count(AADD_MULT_CACHED);
            // 6. multiply w/ product of root weights
            AADD_WGT prod = wgt_mul(AADD_WEIGHT(mat), AADD_WEIGHT(vec));
            AADD_WGT new_weight = wgt_mul(prod, AADD_WEIGHT(res));
            res = aadd_bundle(AADD_TARGET(res), new_weight);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both AADDs
    BDDVAR var;
    AADD vec_low, vec_high, mat_low, mat_high, u00, u10, u01, u11;
    aadd_get_topvar(vec, nextvar, &var, &vec_low, &vec_high);
    aadd_get_topvar(mat, 2*nextvar, &var, &mat_low, &mat_high);
    aadd_get_topvar(mat_low, 2*nextvar+1, &var, &u00, &u10);
    aadd_get_topvar(mat_high,2*nextvar+1, &var, &u01, &u11);

    // 2. propagate "in-between" weights of matrix AADD
    u00 = aadd_bundle(AADD_TARGET(u00), wgt_mul(AADD_WEIGHT(u00), AADD_WEIGHT(mat_low)));
    u10 = aadd_bundle(AADD_TARGET(u10), wgt_mul(AADD_WEIGHT(u10), AADD_WEIGHT(mat_low)));
    u01 = aadd_bundle(AADD_TARGET(u01), wgt_mul(AADD_WEIGHT(u01), AADD_WEIGHT(mat_high)));
    u11 = aadd_bundle(AADD_TARGET(u11), wgt_mul(AADD_WEIGHT(u11), AADD_WEIGHT(mat_high)));

    // 3. recursive calls (4 tasks: SPAWN 3, CALL 1)
    // |u00 u01| |vec_low | = vec_low|u00| + vec_high|u01|
    // |u10 u11| |vec_high|          |u10|           |u11|
    AADD res_low00, res_low10, res_high01, res_high11; //                                       [GC refs stack]
    aadd_refs_spawn(SPAWN(aadd_matvec_mult_rec, u00, vec_low,  nvars, nextvar+1));  // fork 1
    aadd_refs_spawn(SPAWN(aadd_matvec_mult_rec, u10, vec_low,  nvars, nextvar+1));  // fork 2
    aadd_refs_spawn(SPAWN(aadd_matvec_mult_rec, u01, vec_high, nvars, nextvar+1));  // fork 3
    res_high11 = aadd_refs_push(CALL(aadd_matvec_mult_rec, u11, vec_high, nvars, nextvar+1));// [res_high11]
    res_high01 = aadd_refs_sync(SYNC(aadd_matvec_mult_rec));                        // join 3   [res_high11]
    aadd_refs_pop(1);                                                               //          []
    AADD res_high = aadd_refs_push(aadd_makenode(nextvar, res_high01, res_high11)); //          [res_high]
    res_low10  = aadd_refs_push(aadd_refs_sync(SYNC(aadd_matvec_mult_rec)));        // join 2   [res_low10,res_high]
    res_low00  = aadd_refs_sync(SYNC(aadd_matvec_mult_rec));                        // join 1   [res_low10,res_high]
    aadd_refs_pop(1);                                                               //          [res_high]
    AADD res_low  = aadd_refs_push(aadd_makenode(nextvar, res_low00,  res_low10));  //          [res_low,res_high]

    // 4. add resulting AADDs
    res = CALL(aadd_plus, res_low, res_high);                                       //          [res_low,res_high]
    aadd_refs_pop(2);                                                               //          []

    // Insert in cache (before multiplication w/ root weights)
    if (cachenow) {
        if (cache_put3(CACHE_AADD_MATVEC_MULT, nextvar, AADD_TARGET(mat), AADD_TARGET(vec), res)) 
            sylvan_stats_count(AADD_MULT_CACHEDPUT);
    }

    // 5. multiply w/ product of root weights
    AADD_WGT prod = wgt_mul(AADD_WEIGHT(mat), AADD_WEIGHT(vec));
    AADD_WGT new_weight = wgt_mul(prod, AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_weight);

    return res;
}

TASK_IMPL_4(AADD, aadd_matmat_mult_rec, AADD, a, AADD, b, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (AADD_WEIGHT(a) == AADD_ZERO || AADD_WEIGHT(b) == AADD_ZERO)
        return aadd_bundle(AADD_TERMINAL, AADD_ZERO);

    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(AADD_TARGET(a) == AADD_TERMINAL);
        assert(AADD_TARGET(b) == AADD_TERMINAL);
        AADD_WGT prod = wgt_mul(AADD_WEIGHT(a), AADD_WEIGHT(b));
        return aadd_bundle(AADD_TERMINAL, prod);
    }

    sylvan_gc_test();

    // Check cache
    AADD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_AADD_MATMAT_MULT, nextvar, AADD_TARGET(a), AADD_TARGET(b), &res)) {
            sylvan_stats_count(AADD_MULT_CACHED);
            // 7. multiply w/ product of root weights
            AADD_WGT prod = wgt_mul(AADD_WEIGHT(a), AADD_WEIGHT(b));
            AADD_WGT new_weight = wgt_mul(prod, AADD_WEIGHT(res));
            res = aadd_bundle(AADD_TARGET(res), new_weight);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both AADDs
    BDDVAR var;
    AADD a_low, a_high, a00, a10, a01, a11, b_low, b_high, b00, b10, b01, b11;
    aadd_get_topvar(a, 2*nextvar, &var, &a_low, &a_high);
    aadd_get_topvar(b, 2*nextvar, &var, &b_low, &b_high);
    aadd_get_topvar(a_low, 2*nextvar+1, &var, &a00, &a10);
    aadd_get_topvar(a_high,2*nextvar+1, &var, &a01, &a11);
    aadd_get_topvar(b_low, 2*nextvar+1, &var, &b00, &b10);
    aadd_get_topvar(b_high,2*nextvar+1, &var, &b01, &b11);

    // 2. propagate "in-between" weights down
    a00 = aadd_bundle(AADD_TARGET(a00), wgt_mul(AADD_WEIGHT(a_low), AADD_WEIGHT(a00)));
    a10 = aadd_bundle(AADD_TARGET(a10), wgt_mul(AADD_WEIGHT(a_low), AADD_WEIGHT(a10)));
    a01 = aadd_bundle(AADD_TARGET(a01), wgt_mul(AADD_WEIGHT(a_high),AADD_WEIGHT(a01)));
    a11 = aadd_bundle(AADD_TARGET(a11), wgt_mul(AADD_WEIGHT(a_high),AADD_WEIGHT(a11)));
    b00 = aadd_bundle(AADD_TARGET(b00), wgt_mul(AADD_WEIGHT(b_low), AADD_WEIGHT(b00)));
    b10 = aadd_bundle(AADD_TARGET(b10), wgt_mul(AADD_WEIGHT(b_low), AADD_WEIGHT(b10)));
    b01 = aadd_bundle(AADD_TARGET(b01), wgt_mul(AADD_WEIGHT(b_high),AADD_WEIGHT(b01)));
    b11 = aadd_bundle(AADD_TARGET(b11), wgt_mul(AADD_WEIGHT(b_high),AADD_WEIGHT(b11)));

    // 3. recursive calls (8 tasks: SPAWN 7, CALL 1)
    // |a00 a01| |b00 b01| = b00|a00| + b10|a01| , b01|a00| + b11|a01|
    // |a10 a11| |b10 b11|      |a10|      |a11|      |a10|      |a11|
    AADD a00_b00, a00_b01, a10_b00, a10_b01, a01_b10, a01_b11, a11_b10, a11_b11; //         [GC refs stack]
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a00, b00, nvars, nextvar+1)); // fork 1
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a10, b00, nvars, nextvar+1)); // fork 2
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a00, b01, nvars, nextvar+1)); // fork 3
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a10, b01, nvars, nextvar+1)); // fork 4
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a01, b10, nvars, nextvar+1)); // fork 5
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a11, b10, nvars, nextvar+1)); // fork 6
    aadd_refs_spawn(SPAWN(aadd_matmat_mult_rec, a01, b11, nvars, nextvar+1)); // fork 7 
    a11_b11 = aadd_refs_push(CALL(aadd_matmat_mult_rec, a11, b11, nvars, nextvar+1)); //    [a11_b11]
    a01_b11 = aadd_refs_sync(SYNC(aadd_matmat_mult_rec));                     // join 7     [a11_b11]
    aadd_refs_pop(1);                                                         //            []
    AADD rh2 = aadd_refs_push(aadd_makenode(2*nextvar+1, a01_b11, a11_b11));  //            [rh2]
    a11_b10 = aadd_refs_push(aadd_refs_sync(SYNC(aadd_matmat_mult_rec)));     // join 6     [a11_b10,rh2]
    a01_b10 = aadd_refs_sync(SYNC(aadd_matmat_mult_rec));                     // join 5     [a11_b10,rh2]
    aadd_refs_pop(1);                                                         //            [rh2]
    AADD lh2 = aadd_refs_push(aadd_makenode(2*nextvar+1, a01_b10, a11_b10));  //            [lh2,rh2]
    a10_b01 = aadd_refs_push(aadd_refs_sync(SYNC(aadd_matmat_mult_rec)));     // join 4     [b10_b01,lh2,rh2]
    a00_b01 = aadd_refs_sync(SYNC(aadd_matmat_mult_rec));                     // join 3     [b10_b01,lh2,rh2]
    aadd_refs_pop(1);                                                         //            [lh2,rh2]
    AADD rh1 = aadd_refs_push(aadd_makenode(2*nextvar+1, a00_b01, a10_b01));  //            [rh1,lh2,rh2]
    a10_b00 = aadd_refs_push(aadd_refs_sync(SYNC(aadd_matmat_mult_rec)));     // join 2     [b10_b00,rh1,lh2,rh2]
    a00_b00 = aadd_refs_sync(SYNC(aadd_matmat_mult_rec));                     // join 1     [b10_b00,rh1,lh2,rh2]
    aadd_refs_pop(1);                                                         //            [rh1,lh2,rh2]
    AADD lh1 = aadd_refs_push(aadd_makenode(2*nextvar+1, a00_b00, a10_b00));  //            [lh1,rh1,lh2,rh2]

    // 4. add resulting AADDs
    AADD lh, rh;
    aadd_refs_spawn(SPAWN(aadd_plus, lh1, lh2));
    rh = CALL(aadd_plus, rh1, rh2);
    aadd_refs_push(rh);                                                       //            [rh,lh1,rh1,lh2,rh2]
    lh = aadd_refs_sync(SYNC(aadd_plus));                                     //            [rh,lh1,rh1,lh2,rh2]
    aadd_refs_pop(5);                                                         //            []

    // 5. put left and right halves of matix together
    res = aadd_makenode(2*nextvar, lh, rh);

    // Insert in cache
    if (cachenow) {
        if (cache_put3(CACHE_AADD_MATMAT_MULT, nextvar, AADD_TARGET(a), AADD_TARGET(b), res)) 
            sylvan_stats_count(AADD_MULT_CACHEDPUT);
    }

    // 6. multiply w/ product of root weights
    AADD_WGT prod = wgt_mul(AADD_WEIGHT(a), AADD_WEIGHT(b));
    AADD_WGT new_weight = wgt_mul(prod, AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_weight);

    return res;
}

TASK_IMPL_4(AADD_WGT, aadd_inner_product, AADD, a, AADD, b, BDDVAR, nvars, BDDVAR, nextvar)
{
    if (AADD_WEIGHT(a) == AADD_ZERO) return AADD_ZERO;
    if (AADD_WEIGHT(b) == AADD_ZERO) return AADD_ZERO;

    // Terminal case: currently aadd_inner_product doesn't skip variables, 
    // so when both point to terminal, both are scalars.
    // TODO: allow for skipping variables (and multiply res w/ 2^{# skipped})
    // (requires adding some wgt_from_int() function in wgt interface)
    if (nextvar == nvars) {
        return wgt_mul(AADD_WEIGHT(a), wgt_conj(AADD_WEIGHT(b)));
    }

    // Get var(a) and var(b)
    AADD low_a, low_b, high_a, high_b;
    BDDVAR topvar;

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    aadd_get_topvar(a, nextvar, &topvar, &low_a, &high_a);
    aadd_get_topvar(b, nextvar, &topvar, &low_b, &high_b);

    // Check cache
    // TODO: norm cache key? (<a|b> = <b|a>^\dagger)
    AADD_WGT res;
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        if (cache_get4(CACHE_AADD_INPROD, AADD_TARGET(a), AADD_TARGET(b), nextvar, nvars, &res)) {
            res = wgt_mul(res, AADD_WEIGHT(a));
            res = wgt_mul(res, wgt_conj(AADD_WEIGHT(b)));
            return res;
        }
    }
    
    // Recursive calls
    aadd_refs_spawn(SPAWN(aadd_inner_product, high_a, high_b, nvars, nextvar+1));
    AADD_WGT res_low = aadd_refs_push(CALL(aadd_inner_product, low_a, low_b, nvars, nextvar+1));
    AADD_WGT res_high = aadd_refs_sync(SYNC(aadd_inner_product));
    aadd_refs_pop(1);

    res = wgt_add(res_low, res_high);

    // Insert in cache (before multiplication w/ root weights)
    if (cachenow) {
        cache_put4(CACHE_AADD_INPROD, AADD_TARGET(a), AADD_TARGET(b), nextvar, nvars, res);
    }

    // Multiply result with product of weights of a and (conjugate of) b
    // (Note that we can compute the complex conjugate of |b> by taking the 
    // complex conjugate of all edge weights separately, since 
    // (w1 • w2)* = (w2* • w1*) and for scalars (w2* • w1*) = (w1* • w2*).)
    res = wgt_mul(res, AADD_WEIGHT(a));
    res = wgt_mul(res, wgt_conj(AADD_WEIGHT(b)));
    return res;
}

AADD
aadd_increase_all_vars(AADD a, int k)
{
    if (AADD_TARGET(a) == AADD_TERMINAL) {
        return a;
    }

    // Check cache
    AADD res;
    if (cache_get3(CACHE_AADD_INC_VARS, AADD_TARGET(a), k, 0, &res)) {
        return aadd_bundle(res, AADD_WEIGHT(a));
    }

    // Get node info
    AADD low, high;
    aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
    aaddnode_getchilderen(node, &low, &high);
    BDDVAR curvar = aaddnode_getvar(node);

    // Recursive apply to children (TODO: lace?)
    low  = aadd_increase_all_vars(low, k);
    high = aadd_increase_all_vars(high, k);
    res = aadd_makenode(curvar+k, low, high);

    // We assume the input AADDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(AADD_WEIGHT(res) == AADD_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_AADD_INC_VARS, AADD_TARGET(a), k, 0, AADD_TARGET(res));
    return aadd_bundle(AADD_TARGET(res), AADD_WEIGHT(a));
}

AADD
aadd_replace_terminal(AADD a, AADD_TARG b)
{
    if (AADD_TARGET(a) == AADD_TERMINAL) {
        return aadd_bundle(b, AADD_WEIGHT(a));
    }

    // Check cache
    AADD res;
    if (cache_get3(CACHE_AADD_REPLACE_TERMINAL, AADD_TARGET(a), b, 0, &res)) {
        return aadd_bundle(res, AADD_WEIGHT(a));
    }

    // Get node info
    AADD low, high;
    aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
    aaddnode_getchilderen(node, &low, &high);
    BDDVAR var = aaddnode_getvar(node);

    // Recursive apply to children and makenode (TODO: lace?)
    low  = aadd_replace_terminal(low, b);
    high = aadd_replace_terminal(high, b);
    res = aadd_makenode(var, low, high);
    
    // We assume the input AADDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(AADD_WEIGHT(res) == AADD_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_AADD_REPLACE_TERMINAL, AADD_TARGET(a), b, 0, AADD_TARGET(res));
    return aadd_bundle(AADD_TARGET(res), AADD_WEIGHT(a));
}

AADD
aadd_tensor_prod(AADD a, AADD b, BDDVAR nvars_a)
{
    // Shift all vars of 'b' by 'nvars_a'
    b = aadd_increase_all_vars(b, nvars_a);
 
    // Stack 'a' on top of 'b' (and multiply root edge of 'a' with that of 'b')
    AADD res = aadd_replace_terminal(a, AADD_TARGET(b));
    AADD_WGT new_weight = wgt_mul(AADD_WEIGHT(res), AADD_WEIGHT(b));
    res = aadd_bundle(AADD_TARGET(res), new_weight);

    return res;
}

/*************************</Matrix/vector operations>**************************/





/***************************<AADD utility functions>***************************/


/**
 * Helper function for recursive unmarking
 */
static void
aadd_unmark_rec(AADD a)
{
    if (AADD_TARGET(a) == AADD_TERMINAL) return;
    aaddnode_t n = AADD_GETNODE(AADD_TARGET(a));
    if (!aaddnode_getmark(n)) return;
    aaddnode_setmark(n, 0);
    aadd_unmark_rec(aaddnode_getptrlow(n));
    aadd_unmark_rec(aaddnode_getptrhigh(n));
}

/**
 * Counts nodes in the AADD by marking them.
 */
static uint64_t
aadd_nodecount_mark(AADD a)
{
    if (AADD_TARGET(a) == AADD_TERMINAL) return 0; // don't (repeat) count terminal
    aaddnode_t n = AADD_GETNODE(AADD_TARGET(a));
    if (aaddnode_getmark(n)) return 0;
    aaddnode_setmark(n, 1);
    return 1 + aadd_nodecount_mark(aaddnode_getptrlow(n)) + aadd_nodecount_mark(aaddnode_getptrhigh(n));
}

uint64_t
aadd_countnodes(AADD a)
{
    uint64_t res = aadd_nodecount_mark(a) + 1; // (+ 1 for terminal "node")
    aadd_unmark_rec(a);
    return res;
}

/**************************</AADD utility functions>***************************/





/**************************<Printing & file writing>***************************/

/**
 * Pretty prints the information contained in `n`.
 */
static void aaddnode_pprint(aaddnode_t n)
{
    BDDVAR var = aaddnode_getvar(n);
    AADD low, high;
    aaddnode_getchilderen(n, &low, &high);
    printf("[var=%d, low=%" PRIu64 ", high=%" PRIu64 ", ",
             var,
             AADD_TARGET(low),
             AADD_TARGET(high));
    if(AADD_WEIGHT(low) == AADD_ZERO)      printf("a=AADD_ZERO, ");
    else if(AADD_WEIGHT(high) == AADD_ONE)  printf("a=AADD_ONE, ");
    else {
        printf("a=%" PRIu64 ", ", AADD_WEIGHT(low));
        printf("("); wgt_fprint(stdout, AADD_WEIGHT(low)); printf(")");
    }                      
    if(AADD_WEIGHT(high) == AADD_ZERO)     printf("b=AADD_ZERO ");
    else if(AADD_WEIGHT(high) == AADD_ONE) printf("b=AADD_ONE, ");
    else {                     
        printf("b=%" PRIu64, AADD_WEIGHT(high));
        printf("("); wgt_fprint(stdout, AADD_WEIGHT(high)); printf(")");
    }
    printf("]\n");
}

void
_print_aadd(AADD a)
{
    if(AADD_TARGET(a) != AADD_TERMINAL){
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
        if(!aaddnode_getmark(node)){
            aaddnode_setmark(node, 1);
            printf("%" PRIu64 "\t", AADD_TARGET(a));
            aaddnode_pprint(node);
            AADD low, high;
            aaddnode_getchilderen(node, &low, &high);
            _print_aadd(low);
            _print_aadd(high);
        }
    }
}

void
aadd_printnodes(AADD a)
{
    printf("root edge: %" PRIu64 ", %" PRIu64 " = ",AADD_TARGET(a), AADD_WEIGHT(a));
    wgt_fprint(stdout, AADD_WEIGHT(a));
    printf("\n");
    _print_aadd(a);
    aadd_unmark_rec(a);
}

static void
aadd_fprintdot_edge_label(FILE *out, AADD_WGT w)
{
    fprintf(out, ", label=\"");
    if (w == AADD_ONE) {}
    else if (w == AADD_ZERO) { fprintf(out, "0"); }
    else if (w == AADD_MIN_ONE) { fprintf(out, "-1"); }
    else { wgt_fprint(out, w); }
    fprintf(out, "\"");
}

static void
aadd_fprintdot_rec(FILE *out, AADD a, bool draw_zeros)
{
    // terminal node
    if(AADD_TARGET(a) == AADD_TERMINAL) return;

    aaddnode_t n = AADD_GETNODE(AADD_TARGET(a));
    if (aaddnode_getmark(n)) return;
    aaddnode_setmark(n, 1);

    // add this node
    fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n",
            AADD_TARGET(a), aaddnode_getvar(n));

    
    // children of this node
    AADD low, high;
    aaddnode_getchilderen(n, &low, &high);
    aadd_fprintdot_rec(out, low, draw_zeros);
    aadd_fprintdot_rec(out, high, draw_zeros);

    // add edge from this node to each child (unless weight 0)
    if (draw_zeros || AADD_WEIGHT(low) != AADD_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed",
                    AADD_TARGET(a), AADD_TARGET(low));
        aadd_fprintdot_edge_label(out, AADD_WEIGHT(low));
        fprintf(out, "];\n");
    }
    if (draw_zeros || AADD_WEIGHT(high) != AADD_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                    AADD_TARGET(a), AADD_TARGET(high));
        aadd_fprintdot_edge_label(out, AADD_WEIGHT(high));
        fprintf(out, "];\n");
    }
}

void
aadd_fprintdot(FILE *out, AADD a, bool draw_zeros)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir=forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid", AADD_TARGET(a));
    aadd_fprintdot_edge_label(out, AADD_WEIGHT(a));
    fprintf(out, "];\n");

    // terminal node
    fprintf(out, "%" PRIu64 " [shape=box, label=\"T\"];\n", AADD_TERMINAL);

    // recursively add nodes
    aadd_fprintdot_rec(out, a, draw_zeros);
    aadd_unmark_rec(a);

    fprintf(out, "}\n");
}

/**************************</Printing & file writing>**************************/





/********************************<Debug stuff>*********************************/

bool
aadd_equivalent(AADD a, AADD b, int n, bool exact, bool verbose)
{
    bool has_next = true;
    AADD_WGT wgt_a, wgt_b;
    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    while(has_next){
        wgt_a = aadd_getvalue(a, x);
        wgt_b = aadd_getvalue(b, x);
        if(exact){
            if(!wgt_eq(wgt_a, wgt_b)){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", wgt a ="); wgt_fprint(stdout, wgt_a);
                    printf(" != wgt b ="); wgt_fprint(stdout, wgt_b);
                    printf("\n");
                }
                return false;
            }
        }
        else{
            if(!wgt_approx_eq(wgt_a, wgt_b)){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", wgt a ="); wgt_fprint(stdout, wgt_a);
                    printf(" !~= wgt b ="); wgt_fprint(stdout, wgt_b);
                    printf("\n");
                }
                return false;
            }
        }
        has_next = _next_bitstring(x, n);
    }
    return true;
}

bool
aadd_is_ordered_rec(AADD a, BDD nextvar, BDD nvars)
{
    // Terminal case
    if (AADD_TARGET(a) == AADD_TERMINAL) return true;

    // Get top node
    BDDVAR var;
    AADD low, high;
    aadd_get_topvar(a, AADD_INVALID_VAR, &var, &low, &high);

    // Check variable order
    if (var >= nvars || var < nextvar) return false;

    // Check cache
    uint64_t res;
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_AADD_IS_ORDERED, AADD_TARGET(a), 0, 0, &res)) {
            return (bool) res;
        }
    }

    // Recursive calls
    bool res_low, res_high;
    var++;
    res_low = aadd_is_ordered_rec(low, var, nvars);
    res_high = aadd_is_ordered_rec(high, var, nvars);
    var--;
    res = (res_low && res_high);

    // Put res in cache and return
    if (cachenow) cache_put3(CACHE_AADD_IS_ORDERED, AADD_TARGET(a), 0, 0, res);
    return res;
}

bool
_next_bitstring(bool *x, int n)
{
    // binary add 1
    bool success = false;
    for(int k=0; k<n; k++){
        if(x[k] == 1){
            x[k] = 0;
        }
        else if(x[k] == 0) {
            x[k] = 1;
            success = true;
            break;
        }
    }
    return success;
}

void
_print_bitstring(bool *x, int n, bool backwards)
{
    if (backwards)
        for(int k=n-1; k>=0; k--) printf("%d", x[k]);
    else 
        for(int k=0; k<n; k++) printf("%d", x[k]);
}

uint64_t
bitarray_to_int(bool *x, int n, bool MSB_first)
{
    uint64_t res = 0, k = 1;
    if (MSB_first) {
        for (int i = n-1; i >= 0; i--) {
            if (x[i] == 1) res |= k;
            k = k<<1;
        }
    }
    else {
        for (int i = 0; i < n; i++) {
            if (x[i] == 1) res |= k;
            k = k<<1;
        }
    }
    return res;
}

bool *
int_to_bitarray(uint64_t n, int length, bool MSB_first)
{
    if (length <= 0) length = (int) ceil(log2(n));
    bool *res = malloc(sizeof(bool) * length);
    if (MSB_first) {
        for (int i = length-1; i >= 0; i--) {
            res[i] = (n & 1LL);
            n = n>>1;
        }
    }
    else {
        for (int i = 0; i < length; i++) {
            res[i] = (n & 1LL);
            n = n>>1;
        }
    }
    return res;
}

bool
bit_from_int(uint64_t a, uint8_t index)
{
    // assumes index=0 is the LSB
    uint64_t mask = 1<<index;
    uint64_t res = a & mask;
    res = res>>index;
    return (bool) res;
}

void
reverse_bit_array(bool *x, int length)
{
    bool tmp;
    for(int i = 0; i<(length/2); i++){
        tmp = x[i];
        x[i] = x[length-i-1];
        x[length-i-1] = tmp;
    }
}

/*******************************</Debug stuff>*********************************/
