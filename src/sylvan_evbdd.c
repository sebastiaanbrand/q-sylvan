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
#include <sylvan_evbdd.h>
#include <sylvan_refs.h>

static int granularity = 1; // operation cache access granularity


bool larger_wgt_indices;
int weight_norm_strat;
EVBDD_WGT (*normalize_weights)(EVBDD_WGT *, EVBDD_WGT *);

/*******************<Garbage collection, references, marking>******************/

/**
 * Most of this gc code is copy-paste from sylvan_mtbdd.c, however because the 
 * bit structure of BDDs and EVBDDs are a bit different we can't use the mtbdd
 * code directly. Since the way sylvan_mtbdd/sylvan_common/sylvan_table are
 * structured we need to copy-paste a few more functions/variables than we 
 * actually change.
 */

/* 
 * Recursively mark EVBDD nodes as 'in use'.
 * This is really the only gc function which is different for EVBDDs vs MTBDDs.
 */
VOID_TASK_IMPL_1(evbdd_gc_mark_rec, EVBDD, a)
{
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) return;

    if (llmsset_mark(nodes, EVBDD_TARGET(a))) {
        evbddnode_t n = EVBDD_GETNODE(EVBDD_TARGET(a));
        SPAWN(evbdd_gc_mark_rec, evbddnode_getptrlow(n));
        CALL(evbdd_gc_mark_rec, evbddnode_getptrhigh(n));
        SYNC(evbdd_gc_mark_rec);
    }
}

/**
 * External references
 */
refs_table_t evbdd_refs;
refs_table_t evbdd_protected;
static int evbdd_protected_created = 0;

void
evbdd_protect(EVBDD *a)
{
    if (!evbdd_protected_created) {
        // In C++, sometimes mtbdd_protect is called before Sylvan is initialized. Just create a table.
        protect_create(&evbdd_protected, 4096);
        evbdd_protected_created = 1;
    }
    protect_up(&evbdd_protected, (size_t)a);
}

void
evbdd_unprotect(EVBDD *a)
{
    if (evbdd_protected.refs_table != NULL) protect_down(&evbdd_protected, (size_t)a);
}

size_t
evbdd_count_protected()
{
    return protect_count(&evbdd_protected);
}

/* Called during garbage collection */
VOID_TASK_0(evbdd_gc_mark_external_refs)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = refs_iter(&evbdd_refs, 0, evbdd_refs.refs_size);
    while (it != NULL) {
        SPAWN(evbdd_gc_mark_rec, refs_next(&evbdd_refs, &it, evbdd_refs.refs_size));
        count++;
    }
    while (count--) {
        SYNC(evbdd_gc_mark_rec);
    }
}

/* Called during garbage collection */
VOID_TASK_0(evbdd_gc_mark_protected)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = protect_iter(&evbdd_protected, 0, evbdd_protected.refs_size);
    while (it != NULL) {
        EVBDD *to_mark = (EVBDD*)protect_next(&evbdd_protected, &it, evbdd_protected.refs_size);
        SPAWN(evbdd_gc_mark_rec, *to_mark);
        count++;
    }
    while (count--) {
        SYNC(evbdd_gc_mark_rec);
    }
}

/* Infrastructure for internal markings */
typedef struct evbdd_refs_task
{
    Task *t;
    void *f;
} *evbdd_refs_task_t;

typedef struct evbdd_refs_internal
{
    const EVBDD **pbegin, **pend, **pcur;
    EVBDD *rbegin, *rend, *rcur;
    evbdd_refs_task_t sbegin, send, scur;
} *evbdd_refs_internal_t;

DECLARE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);

VOID_TASK_2(evbdd_refs_mark_p_par, const EVBDD**, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            evbdd_gc_mark_rec(**(begin++));
            count--;
        }
    } else {
        SPAWN(evbdd_refs_mark_p_par, begin, count / 2);
        CALL(evbdd_refs_mark_p_par, begin + (count / 2), count - count / 2);
        SYNC(evbdd_refs_mark_p_par);
    }
}

VOID_TASK_2(evbdd_refs_mark_r_par, EVBDD*, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            evbdd_gc_mark_rec(*begin++);
            count--;
        }
    } else {
        SPAWN(evbdd_refs_mark_r_par, begin, count / 2);
        CALL(evbdd_refs_mark_r_par, begin + (count / 2), count - count / 2);
        SYNC(evbdd_refs_mark_r_par);
    }
}

VOID_TASK_2(evbdd_refs_mark_s_par, evbdd_refs_task_t, begin, size_t, count)
{
    if (count < 32) {
        while (count > 0) {
            Task *t = begin->t;
            if (!TASK_IS_STOLEN(t)) return;
            if (t->f == begin->f && TASK_IS_COMPLETED(t)) {
                evbdd_gc_mark_rec(*(EVBDD*)TASK_RESULT(t));
            }
            begin += 1;
            count -= 1;
        }
    } else {
        if (!TASK_IS_STOLEN(begin->t)) return;
        SPAWN(evbdd_refs_mark_s_par, begin, count / 2);
        CALL(evbdd_refs_mark_s_par, begin + (count / 2), count - count / 2);
        SYNC(evbdd_refs_mark_s_par);
    }
}

VOID_TASK_0(evbdd_refs_mark_task)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    SPAWN(evbdd_refs_mark_p_par, evbdd_refs_key->pbegin, evbdd_refs_key->pcur-evbdd_refs_key->pbegin);
    SPAWN(evbdd_refs_mark_r_par, evbdd_refs_key->rbegin, evbdd_refs_key->rcur-evbdd_refs_key->rbegin);
    CALL(evbdd_refs_mark_s_par, evbdd_refs_key->sbegin, evbdd_refs_key->scur-evbdd_refs_key->sbegin);
    SYNC(evbdd_refs_mark_r_par);
    SYNC(evbdd_refs_mark_p_par);
}

/* Called during garbage collection */
VOID_TASK_0(evbdd_refs_mark)
{
    TOGETHER(evbdd_refs_mark_task);
}

VOID_TASK_0(evbdd_refs_init_task)
{
    evbdd_refs_internal_t s = (evbdd_refs_internal_t)malloc(sizeof(struct evbdd_refs_internal));
    s->pcur = s->pbegin = (const EVBDD**)malloc(sizeof(EVBDD*) * 1024);
    s->pend = s->pbegin + 1024;
    s->rcur = s->rbegin = (EVBDD*)malloc(sizeof(EVBDD) * 1024);
    s->rend = s->rbegin + 1024;
    s->scur = s->sbegin = (evbdd_refs_task_t)malloc(sizeof(struct evbdd_refs_task) * 1024);
    s->send = s->sbegin + 1024;
    SET_THREAD_LOCAL(evbdd_refs_key, s);
}

VOID_TASK_0(evbdd_refs_init)
{
    INIT_THREAD_LOCAL(evbdd_refs_key);
    TOGETHER(evbdd_refs_init_task);
    sylvan_gc_add_mark(TASK(evbdd_refs_mark));
}

VOID_TASK_0(evbdd_refs_cleanup_task)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    free(evbdd_refs_key->pbegin);
    free(evbdd_refs_key->rbegin);
    free(evbdd_refs_key->sbegin);
    free(evbdd_refs_key);
}

/**
 * Called by evbdd_quit. Cleans up the (thread local) malloc'ed evbdd_refs_key.
 * 
 * NOTE: this cleanup isn't done in sylvan_mtbdd.c, but not doing this causes
 * memory leaks when calling initializing and quiting Sylvan multiple times
 * during the same program run.
 */
VOID_TASK_0(evbdd_refs_cleanup)
{
    TOGETHER(evbdd_refs_cleanup_task);
}

void
evbdd_refs_ptrs_up(evbdd_refs_internal_t evbdd_refs_key)
{
    size_t cur = evbdd_refs_key->pcur - evbdd_refs_key->pbegin;
    size_t size = evbdd_refs_key->pend - evbdd_refs_key->pbegin;
    evbdd_refs_key->pbegin = (const EVBDD**)realloc(evbdd_refs_key->pbegin, sizeof(EVBDD*) * size * 2);
    evbdd_refs_key->pcur = evbdd_refs_key->pbegin + cur;
    evbdd_refs_key->pend = evbdd_refs_key->pbegin + (size * 2);
}

EVBDD __attribute__((noinline))
evbdd_refs_refs_up(evbdd_refs_internal_t evbdd_refs_key, EVBDD res)
{
    long size = evbdd_refs_key->rend - evbdd_refs_key->rbegin;
    evbdd_refs_key->rbegin = (EVBDD*)realloc(evbdd_refs_key->rbegin, sizeof(EVBDD) * size * 2);
    evbdd_refs_key->rcur = evbdd_refs_key->rbegin + size;
    evbdd_refs_key->rend = evbdd_refs_key->rbegin + (size * 2);
    return res;
}

void __attribute__((noinline))
evbdd_refs_tasks_up(evbdd_refs_internal_t evbdd_refs_key)
{
    long size = evbdd_refs_key->send - evbdd_refs_key->sbegin;
    evbdd_refs_key->sbegin = (evbdd_refs_task_t)realloc(evbdd_refs_key->sbegin, sizeof(struct evbdd_refs_task) * size * 2);
    evbdd_refs_key->scur = evbdd_refs_key->sbegin + size;
    evbdd_refs_key->send = evbdd_refs_key->sbegin + (size * 2);
}

void __attribute__((unused))
evbdd_refs_pushptr(const EVBDD *ptr)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    *evbdd_refs_key->pcur++ = ptr;
    if (evbdd_refs_key->pcur == evbdd_refs_key->pend) evbdd_refs_ptrs_up(evbdd_refs_key);
}

void __attribute__((unused))
evbdd_refs_popptr(size_t amount)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    evbdd_refs_key->pcur -= amount;
}

EVBDD __attribute__((unused))
evbdd_refs_push(EVBDD a)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    *(evbdd_refs_key->rcur++) = a;
    if (evbdd_refs_key->rcur == evbdd_refs_key->rend) return evbdd_refs_refs_up(evbdd_refs_key, a);
    else return a;
}

void __attribute__((unused))
evbdd_refs_pop(long amount)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    evbdd_refs_key->rcur -= amount;
}

void
evbdd_refs_spawn(Task *t)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    evbdd_refs_key->scur->t = t;
    evbdd_refs_key->scur->f = t->f;
    evbdd_refs_key->scur += 1;
    if (evbdd_refs_key->scur == evbdd_refs_key->send) evbdd_refs_tasks_up(evbdd_refs_key);
}

EVBDD
evbdd_refs_sync(EVBDD result)
{
    LOCALIZE_THREAD_LOCAL(evbdd_refs_key, evbdd_refs_internal_t);
    evbdd_refs_key->scur -= 1;
    return result;
}

/******************</Garbage collection, references, marking>******************/





/*************************<Cleaning edge weight table>*************************/

static int auto_gc_wgt_table  = 1;
static double wgt_table_gc_thres = 0.5;

void
evbdd_set_auto_gc_wgt_table(bool enabled)
{
    auto_gc_wgt_table = enabled;
}

void
evbdd_set_gc_wgt_table_thres(double fraction_filled)
{
    wgt_table_gc_thres = fraction_filled;
}

double
evbdd_get_gc_wgt_table_thres()
{
    return wgt_table_gc_thres;
}


void
evbdd_gc_wgt_table()
{
    // gc edge weight table and keep wgts of protected EVBDDs (and update those)
    // 1. Create new edge weight table table
    wgt_table_gc_init_new(init_wgt_table_entries);

    // 2. Fill new table with wgts in protected EVBDDs and update those EVBDDs
    uint64_t *it = protect_iter(&evbdd_protected, 0, evbdd_protected.refs_size);
    while (it != NULL) {
        EVBDD *to_protect_wgts = (EVBDD*)protect_next(&evbdd_protected, &it, evbdd_protected.refs_size);
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

TASK_IMPL_1(EVBDD, _fill_new_wgt_table, EVBDD, a)
{
    // Check cache
    EVBDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_EVBDD_CLEAN_WGT_TABLE, 0LL, a, 0LL, &res)) {
            return res;
        }
    }

    // Move weight from old to new table, get new index
    EVBDD_WGT new_wgt = wgt_table_gc_keep(EVBDD_WEIGHT(a));
    a = evbdd_bundle(EVBDD_TARGET(a), new_wgt);

    // If terminal, return
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) return a;
    
    // Recursive for children
    EVBDD low, high;
    evbddnode_t n = EVBDD_GETNODE(EVBDD_TARGET(a));
    evbddnode_getchilderen(n, &low, &high);
    evbdd_refs_spawn(SPAWN(_fill_new_wgt_table, high));
    low = CALL(_fill_new_wgt_table, low);
    evbdd_refs_push(low);
    high = evbdd_refs_sync(SYNC(_fill_new_wgt_table));
    evbdd_refs_pop(1);

    // We don't need to use the 'evbdd_makenode()' function which normalizes the 
    // weights, because the EVBDD doesn't actually change, only the WGT indices,
    // but none of the actual values.
    EVBDD_TARG ptr = _evbdd_makenode(evbddnode_getvar(n), EVBDD_TARGET(low), EVBDD_TARGET(high), EVBDD_WEIGHT(low), EVBDD_WEIGHT(high));

    // Put in cache, return
    res = evbdd_bundle(ptr, new_wgt);
    if (cachenow) cache_put3(CACHE_EVBDD_CLEAN_WGT_TABLE, 0LL, a, 0LL, res);
    return res;
}

bool
evbdd_test_gc_wgt_table()
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
static int evbdd_initialized = 0;

static void
evbdd_quit()
{
    refs_free(&evbdd_refs);
    if (evbdd_protected_created) {
        protect_free(&evbdd_protected);
        evbdd_protected_created = 0;
    }
    RUN(evbdd_refs_cleanup);
    evbdd_initialized = 0;
    sylvan_edge_weights_free();
}

void
sylvan_init_evbdd(size_t min_wgt_tablesize, size_t max_wgt_tablesize,
                 double wgt_tab_tolerance, int edge_weigth_backend, 
                 int norm_strat, void *init_wgt_tab_entries)
{
    if (evbdd_initialized) return;
    evbdd_initialized = 1;

    int index_size = (int) ceil(log2(max_wgt_tablesize));
    if (index_size > 33) {
        fprintf(stderr,"max edge weight storage size is 2^33 (2^16 when using storing r and i seperately)\n");
        exit(1);
    }
    if (index_size > 23) larger_wgt_indices = true;
    else larger_wgt_indices = false;

    sylvan_register_quit(evbdd_quit);
    sylvan_gc_add_mark(TASK(evbdd_gc_mark_external_refs));
    sylvan_gc_add_mark(TASK(evbdd_gc_mark_protected));

    refs_create(&evbdd_refs, 1024);
    if (!evbdd_protected_created) {
        protect_create(&evbdd_protected, 4096);
        evbdd_protected_created = 1;
    }

    // TODO: pass edge weight type to sylvan_init_evbdd
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

    RUN(evbdd_refs_init);
}

void
sylvan_init_evbdd_defaults(size_t min_wgt_tablesize, size_t max_wgt_tablesize)
{
    sylvan_init_evbdd(min_wgt_tablesize, max_wgt_tablesize, -1, COMP_HASHMAP, NORM_LOW, NULL);
}


void
evbdd_set_caching_granularity(int g)
{
    granularity = g;
}

/*****************************</Initialization>********************************/





/**************************<Matrix/vector operations>**************************/

EVBDD_WGT
evbdd_getvalue(EVBDD a, bool* path)
{
    EVBDD_WGT res = EVBDD_ONE;
    EVBDD low, high;
    for (;;) {
        res = wgt_mul(res, EVBDD_WEIGHT(a));
        
        // if the current edge is pointing to the terminal, we're done.
        if (EVBDD_TARGET(a) == EVBDD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
        BDDVAR var     = evbddnode_getvar(node);
        evbddnode_getchilderen(node, &low, &high);

        // Condition low/high choice on basis state vector[var]
        a = (path[var] == 0) ? low : high;
    }

    return res;
}

static void
evbdd_do_before_mult(EVBDD *a, EVBDD *b)
{
    // check if edge weight table needs gc
    if (auto_gc_wgt_table && evbdd_test_gc_wgt_table()) {
        evbdd_protect(a);
        evbdd_protect(b);
        evbdd_gc_wgt_table();
        evbdd_unprotect(a);
        evbdd_unprotect(b);
    }
}

static void
norm_commuting_cache_key(EVBDD a, EVBDD b, EVBDD *x, EVBDD *y)
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

TASK_IMPL_2(EVBDD, evbdd_plus, EVBDD, a, EVBDD, b)
{
    // Trivial cases
    if(EVBDD_WEIGHT(a) == EVBDD_ZERO) return b;
    if(EVBDD_WEIGHT(b) == EVBDD_ZERO) return a;

    sylvan_gc_test();

    // Get var(a) and var(b)
    EVBDD low_a, low_b, high_a, high_b, res;
    BDDVAR var_a = UINT32_MAX, var_b = UINT32_MAX, topvar;
    if (EVBDD_TARGET(a) != EVBDD_TERMINAL) {
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
        var_a  = evbddnode_getvar(node);
    }
    if (EVBDD_TARGET(b) != EVBDD_TERMINAL) {
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(b));
        var_b  = evbddnode_getvar(node);
    }

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    evbdd_get_topvar(a, var_b, &topvar, &low_a, &high_a);
    evbdd_get_topvar(b, var_a, &topvar, &low_b, &high_b);

    // Base/terminal case: same target and same variable
    if(EVBDD_TARGET(a) == EVBDD_TARGET(b) && var_a == var_b){
        EVBDD_WGT sum = wgt_add(EVBDD_WEIGHT(a), EVBDD_WEIGHT(b));
        res = evbdd_bundle(EVBDD_TARGET(a), sum);
        return res;
    }

    // Check cache
    EVBDD x, y;
    norm_commuting_cache_key(a, b, &x, &y); // (a+b) = (b+a) so normalize cache key
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_EVBDD_PLUS, sylvan_false, x, y, &res)) {
            sylvan_stats_count(EVBDD_PLUS_CACHED);
            return res;
        }
    }

    // If not base/terminal case, pass edge weight of current edge down
    EVBDD_WGT wgt_la, wgt_ha, wgt_lb, wgt_hb;
    wgt_la = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(low_a));
    wgt_ha = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(high_a));
    wgt_lb = wgt_mul(EVBDD_WEIGHT(b), EVBDD_WEIGHT(low_b));
    wgt_hb = wgt_mul(EVBDD_WEIGHT(b), EVBDD_WEIGHT(high_b));
    low_a  = evbdd_refs_push(evbdd_bundle(EVBDD_TARGET(low_a),  wgt_la));
    high_a = evbdd_refs_push(evbdd_bundle(EVBDD_TARGET(high_a), wgt_ha));
    low_b  = evbdd_refs_push(evbdd_bundle(EVBDD_TARGET(low_b),  wgt_lb));
    high_b = evbdd_refs_push(evbdd_bundle(EVBDD_TARGET(high_b), wgt_hb));

    // Recursive calls down
    evbdd_refs_spawn(SPAWN(evbdd_plus, high_a, high_b));
    EVBDD low = evbdd_refs_push(CALL(evbdd_plus, low_a, low_b));
    EVBDD high = evbdd_refs_sync(SYNC(evbdd_plus));
    evbdd_refs_pop(5);

    // Put in cache, return
    res = evbdd_makenode(topvar, low, high);
    if (cachenow) {
        if (cache_put3(CACHE_EVBDD_PLUS, sylvan_false, x, y, res)) 
            sylvan_stats_count(EVBDD_PLUS_CACHEDPUT);
    }
    return res;
}



/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(EVBDD, evbdd_matvec_mult, EVBDD, mat, EVBDD, vec, BDDVAR, nvars)
{
    evbdd_do_before_mult(&mat, &vec);
    evbdd_refs_push(mat); evbdd_refs_push(vec);
    EVBDD res = CALL(evbdd_matvec_mult_rec, mat, vec, nvars, 0);
    evbdd_refs_pop(2);
    return res;
}

/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(EVBDD, evbdd_matmat_mult, EVBDD, a, EVBDD, b, BDDVAR, nvars)
{
    evbdd_do_before_mult(&a, &b);
    evbdd_refs_push(a); evbdd_refs_push(b);
    EVBDD res = CALL(evbdd_matmat_mult_rec, a, b, nvars, 0);
    evbdd_refs_pop(2);
    return res;
}

TASK_IMPL_4(EVBDD, evbdd_matvec_mult_rec, EVBDD, mat, EVBDD, vec, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (EVBDD_WEIGHT(mat) == EVBDD_ZERO || EVBDD_WEIGHT(vec) == EVBDD_ZERO)
        return evbdd_bundle(EVBDD_TERMINAL, EVBDD_ZERO);
    
    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(EVBDD_TARGET(mat) == EVBDD_TERMINAL);
        assert(EVBDD_TARGET(vec) == EVBDD_TERMINAL);
        EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(mat), EVBDD_WEIGHT(vec));
        return evbdd_bundle(EVBDD_TERMINAL, prod);
    }

    sylvan_gc_test();

    // Check cache
    EVBDD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_EVBDD_MATVEC_MULT, nextvar, EVBDD_TARGET(mat), EVBDD_TARGET(vec), &res)) {
            sylvan_stats_count(EVBDD_MULT_CACHED);
            // 6. multiply w/ product of root weights
            EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(mat), EVBDD_WEIGHT(vec));
            EVBDD_WGT new_weight = wgt_mul(prod, EVBDD_WEIGHT(res));
            res = evbdd_bundle(EVBDD_TARGET(res), new_weight);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both EVBDDs
    BDDVAR var;
    EVBDD vec_low, vec_high, mat_low, mat_high, u00, u10, u01, u11;
    evbdd_get_topvar(vec, nextvar, &var, &vec_low, &vec_high);
    evbdd_get_topvar(mat, 2*nextvar, &var, &mat_low, &mat_high);
    evbdd_get_topvar(mat_low, 2*nextvar+1, &var, &u00, &u10);
    evbdd_get_topvar(mat_high,2*nextvar+1, &var, &u01, &u11);

    // 2. propagate "in-between" weights of matrix EVBDD
    u00 = evbdd_bundle(EVBDD_TARGET(u00), wgt_mul(EVBDD_WEIGHT(u00), EVBDD_WEIGHT(mat_low)));
    u10 = evbdd_bundle(EVBDD_TARGET(u10), wgt_mul(EVBDD_WEIGHT(u10), EVBDD_WEIGHT(mat_low)));
    u01 = evbdd_bundle(EVBDD_TARGET(u01), wgt_mul(EVBDD_WEIGHT(u01), EVBDD_WEIGHT(mat_high)));
    u11 = evbdd_bundle(EVBDD_TARGET(u11), wgt_mul(EVBDD_WEIGHT(u11), EVBDD_WEIGHT(mat_high)));

    // 3. recursive calls (4 tasks: SPAWN 3, CALL 1)
    // |u00 u01| |vec_low | = vec_low|u00| + vec_high|u01|
    // |u10 u11| |vec_high|          |u10|           |u11|
    EVBDD res_low00, res_low10, res_high01, res_high11; //                                       [GC refs stack]
    evbdd_refs_spawn(SPAWN(evbdd_matvec_mult_rec, u00, vec_low,  nvars, nextvar+1));  // fork 1
    evbdd_refs_spawn(SPAWN(evbdd_matvec_mult_rec, u10, vec_low,  nvars, nextvar+1));  // fork 2
    evbdd_refs_spawn(SPAWN(evbdd_matvec_mult_rec, u01, vec_high, nvars, nextvar+1));  // fork 3
    res_high11 = evbdd_refs_push(CALL(evbdd_matvec_mult_rec, u11, vec_high, nvars, nextvar+1));// [res_high11]
    res_high01 = evbdd_refs_sync(SYNC(evbdd_matvec_mult_rec));                        // join 3   [res_high11]
    evbdd_refs_pop(1);                                                               //          []
    EVBDD res_high = evbdd_refs_push(evbdd_makenode(nextvar, res_high01, res_high11)); //          [res_high]
    res_low10  = evbdd_refs_push(evbdd_refs_sync(SYNC(evbdd_matvec_mult_rec)));        // join 2   [res_low10,res_high]
    res_low00  = evbdd_refs_sync(SYNC(evbdd_matvec_mult_rec));                        // join 1   [res_low10,res_high]
    evbdd_refs_pop(1);                                                               //          [res_high]
    EVBDD res_low  = evbdd_refs_push(evbdd_makenode(nextvar, res_low00,  res_low10));  //          [res_low,res_high]

    // 4. add resulting EVBDDs
    res = CALL(evbdd_plus, res_low, res_high);                                       //          [res_low,res_high]
    evbdd_refs_pop(2);                                                               //          []

    // Insert in cache (before multiplication w/ root weights)
    if (cachenow) {
        if (cache_put3(CACHE_EVBDD_MATVEC_MULT, nextvar, EVBDD_TARGET(mat), EVBDD_TARGET(vec), res)) 
            sylvan_stats_count(EVBDD_MULT_CACHEDPUT);
    }

    // 5. multiply w/ product of root weights
    EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(mat), EVBDD_WEIGHT(vec));
    EVBDD_WGT new_weight = wgt_mul(prod, EVBDD_WEIGHT(res));
    res = evbdd_bundle(EVBDD_TARGET(res), new_weight);

    return res;
}

TASK_IMPL_4(EVBDD, evbdd_matmat_mult_rec, EVBDD, a, EVBDD, b, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (EVBDD_WEIGHT(a) == EVBDD_ZERO || EVBDD_WEIGHT(b) == EVBDD_ZERO)
        return evbdd_bundle(EVBDD_TERMINAL, EVBDD_ZERO);

    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(EVBDD_TARGET(a) == EVBDD_TERMINAL);
        assert(EVBDD_TARGET(b) == EVBDD_TERMINAL);
        EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(b));
        return evbdd_bundle(EVBDD_TERMINAL, prod);
    }

    sylvan_gc_test();

    // Check cache
    EVBDD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_EVBDD_MATMAT_MULT, nextvar, EVBDD_TARGET(a), EVBDD_TARGET(b), &res)) {
            sylvan_stats_count(EVBDD_MULT_CACHED);
            // 7. multiply w/ product of root weights
            EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(b));
            EVBDD_WGT new_weight = wgt_mul(prod, EVBDD_WEIGHT(res));
            res = evbdd_bundle(EVBDD_TARGET(res), new_weight);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both EVBDDs
    BDDVAR var;
    EVBDD a_low, a_high, a00, a10, a01, a11, b_low, b_high, b00, b10, b01, b11;
    evbdd_get_topvar(a, 2*nextvar, &var, &a_low, &a_high);
    evbdd_get_topvar(b, 2*nextvar, &var, &b_low, &b_high);
    evbdd_get_topvar(a_low, 2*nextvar+1, &var, &a00, &a10);
    evbdd_get_topvar(a_high,2*nextvar+1, &var, &a01, &a11);
    evbdd_get_topvar(b_low, 2*nextvar+1, &var, &b00, &b10);
    evbdd_get_topvar(b_high,2*nextvar+1, &var, &b01, &b11);

    // 2. propagate "in-between" weights down
    a00 = evbdd_bundle(EVBDD_TARGET(a00), wgt_mul(EVBDD_WEIGHT(a_low), EVBDD_WEIGHT(a00)));
    a10 = evbdd_bundle(EVBDD_TARGET(a10), wgt_mul(EVBDD_WEIGHT(a_low), EVBDD_WEIGHT(a10)));
    a01 = evbdd_bundle(EVBDD_TARGET(a01), wgt_mul(EVBDD_WEIGHT(a_high),EVBDD_WEIGHT(a01)));
    a11 = evbdd_bundle(EVBDD_TARGET(a11), wgt_mul(EVBDD_WEIGHT(a_high),EVBDD_WEIGHT(a11)));
    b00 = evbdd_bundle(EVBDD_TARGET(b00), wgt_mul(EVBDD_WEIGHT(b_low), EVBDD_WEIGHT(b00)));
    b10 = evbdd_bundle(EVBDD_TARGET(b10), wgt_mul(EVBDD_WEIGHT(b_low), EVBDD_WEIGHT(b10)));
    b01 = evbdd_bundle(EVBDD_TARGET(b01), wgt_mul(EVBDD_WEIGHT(b_high),EVBDD_WEIGHT(b01)));
    b11 = evbdd_bundle(EVBDD_TARGET(b11), wgt_mul(EVBDD_WEIGHT(b_high),EVBDD_WEIGHT(b11)));

    // 3. recursive calls (8 tasks: SPAWN 7, CALL 1)
    // |a00 a01| |b00 b01| = b00|a00| + b10|a01| , b01|a00| + b11|a01|
    // |a10 a11| |b10 b11|      |a10|      |a11|      |a10|      |a11|
    EVBDD a00_b00, a00_b01, a10_b00, a10_b01, a01_b10, a01_b11, a11_b10, a11_b11; //         [GC refs stack]
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a00, b00, nvars, nextvar+1)); // fork 1
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a10, b00, nvars, nextvar+1)); // fork 2
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a00, b01, nvars, nextvar+1)); // fork 3
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a10, b01, nvars, nextvar+1)); // fork 4
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a01, b10, nvars, nextvar+1)); // fork 5
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a11, b10, nvars, nextvar+1)); // fork 6
    evbdd_refs_spawn(SPAWN(evbdd_matmat_mult_rec, a01, b11, nvars, nextvar+1)); // fork 7 
    a11_b11 = evbdd_refs_push(CALL(evbdd_matmat_mult_rec, a11, b11, nvars, nextvar+1)); //    [a11_b11]
    a01_b11 = evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec));                     // join 7     [a11_b11]
    evbdd_refs_pop(1);                                                         //            []
    EVBDD rh2 = evbdd_refs_push(evbdd_makenode(2*nextvar+1, a01_b11, a11_b11));  //            [rh2]
    a11_b10 = evbdd_refs_push(evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec)));     // join 6     [a11_b10,rh2]
    a01_b10 = evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec));                     // join 5     [a11_b10,rh2]
    evbdd_refs_pop(1);                                                         //            [rh2]
    EVBDD lh2 = evbdd_refs_push(evbdd_makenode(2*nextvar+1, a01_b10, a11_b10));  //            [lh2,rh2]
    a10_b01 = evbdd_refs_push(evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec)));     // join 4     [b10_b01,lh2,rh2]
    a00_b01 = evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec));                     // join 3     [b10_b01,lh2,rh2]
    evbdd_refs_pop(1);                                                         //            [lh2,rh2]
    EVBDD rh1 = evbdd_refs_push(evbdd_makenode(2*nextvar+1, a00_b01, a10_b01));  //            [rh1,lh2,rh2]
    a10_b00 = evbdd_refs_push(evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec)));     // join 2     [b10_b00,rh1,lh2,rh2]
    a00_b00 = evbdd_refs_sync(SYNC(evbdd_matmat_mult_rec));                     // join 1     [b10_b00,rh1,lh2,rh2]
    evbdd_refs_pop(1);                                                         //            [rh1,lh2,rh2]
    EVBDD lh1 = evbdd_refs_push(evbdd_makenode(2*nextvar+1, a00_b00, a10_b00));  //            [lh1,rh1,lh2,rh2]

    // 4. add resulting EVBDDs
    EVBDD lh, rh;
    evbdd_refs_spawn(SPAWN(evbdd_plus, lh1, lh2));
    rh = CALL(evbdd_plus, rh1, rh2);
    evbdd_refs_push(rh);                                                       //            [rh,lh1,rh1,lh2,rh2]
    lh = evbdd_refs_sync(SYNC(evbdd_plus));                                     //            [rh,lh1,rh1,lh2,rh2]
    evbdd_refs_pop(5);                                                         //            []

    // 5. put left and right halves of matix together
    res = evbdd_makenode(2*nextvar, lh, rh);

    // Insert in cache
    if (cachenow) {
        if (cache_put3(CACHE_EVBDD_MATMAT_MULT, nextvar, EVBDD_TARGET(a), EVBDD_TARGET(b), res)) 
            sylvan_stats_count(EVBDD_MULT_CACHEDPUT);
    }

    // 6. multiply w/ product of root weights
    EVBDD_WGT prod = wgt_mul(EVBDD_WEIGHT(a), EVBDD_WEIGHT(b));
    EVBDD_WGT new_weight = wgt_mul(prod, EVBDD_WEIGHT(res));
    res = evbdd_bundle(EVBDD_TARGET(res), new_weight);

    return res;
}

TASK_IMPL_4(EVBDD_WGT, evbdd_inner_product, EVBDD, a, EVBDD, b, BDDVAR, nvars, BDDVAR, nextvar)
{
    if (EVBDD_WEIGHT(a) == EVBDD_ZERO) return EVBDD_ZERO;
    if (EVBDD_WEIGHT(b) == EVBDD_ZERO) return EVBDD_ZERO;

    // Terminal case: currently evbdd_inner_product doesn't skip variables, 
    // so when both point to terminal, both are scalars.
    // TODO: allow for skipping variables (and multiply res w/ 2^{# skipped})
    // (requires adding some wgt_from_int() function in wgt interface)
    if (nextvar == nvars) {
        return wgt_mul(EVBDD_WEIGHT(a), wgt_conj(EVBDD_WEIGHT(b)));
    }

    // Get var(a) and var(b)
    EVBDD low_a, low_b, high_a, high_b;
    BDDVAR topvar;

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    evbdd_get_topvar(a, nextvar, &topvar, &low_a, &high_a);
    evbdd_get_topvar(b, nextvar, &topvar, &low_b, &high_b);

    // Check cache
    // TODO: norm cache key? (<a|b> = <b|a>^\dagger)
    EVBDD_WGT res;
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        if (cache_get4(CACHE_EVBDD_INPROD, EVBDD_TARGET(a), EVBDD_TARGET(b), nextvar, nvars, &res)) {
            res = wgt_mul(res, EVBDD_WEIGHT(a));
            res = wgt_mul(res, wgt_conj(EVBDD_WEIGHT(b)));
            return res;
        }
    }
    
    // Recursive calls
    evbdd_refs_spawn(SPAWN(evbdd_inner_product, high_a, high_b, nvars, nextvar+1));
    EVBDD_WGT res_low = evbdd_refs_push(CALL(evbdd_inner_product, low_a, low_b, nvars, nextvar+1));
    EVBDD_WGT res_high = evbdd_refs_sync(SYNC(evbdd_inner_product));
    evbdd_refs_pop(1);

    res = wgt_add(res_low, res_high);

    // Insert in cache (before multiplication w/ root weights)
    if (cachenow) {
        cache_put4(CACHE_EVBDD_INPROD, EVBDD_TARGET(a), EVBDD_TARGET(b), nextvar, nvars, res);
    }

    // Multiply result with product of weights of a and (conjugate of) b
    // (Note that we can compute the complex conjugate of |b> by taking the 
    // complex conjugate of all edge weights separately, since 
    // (w1 • w2)* = (w2* • w1*) and for scalars (w2* • w1*) = (w1* • w2*).)
    res = wgt_mul(res, EVBDD_WEIGHT(a));
    res = wgt_mul(res, wgt_conj(EVBDD_WEIGHT(b)));
    return res;
}

EVBDD
evbdd_increase_all_vars(EVBDD a, int k)
{
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) {
        return a;
    }

    // Check cache
    EVBDD res;
    if (cache_get3(CACHE_EVBDD_INC_VARS, EVBDD_TARGET(a), k, 0, &res)) {
        return evbdd_bundle(res, EVBDD_WEIGHT(a));
    }

    // Get node info
    EVBDD low, high;
    evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
    evbddnode_getchilderen(node, &low, &high);
    BDDVAR curvar = evbddnode_getvar(node);

    // Recursive apply to children (TODO: lace?)
    low  = evbdd_increase_all_vars(low, k);
    high = evbdd_increase_all_vars(high, k);
    res = evbdd_makenode(curvar+k, low, high);

    // We assume the input EVBDDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(EVBDD_WEIGHT(res) == EVBDD_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_EVBDD_INC_VARS, EVBDD_TARGET(a), k, 0, EVBDD_TARGET(res));
    return evbdd_bundle(EVBDD_TARGET(res), EVBDD_WEIGHT(a));
}

EVBDD
evbdd_replace_terminal(EVBDD a, EVBDD_TARG b)
{
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) {
        return evbdd_bundle(b, EVBDD_WEIGHT(a));
    }

    // Check cache
    EVBDD res;
    if (cache_get3(CACHE_EVBDD_REPLACE_TERMINAL, EVBDD_TARGET(a), b, 0, &res)) {
        return evbdd_bundle(res, EVBDD_WEIGHT(a));
    }

    // Get node info
    EVBDD low, high;
    evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
    evbddnode_getchilderen(node, &low, &high);
    BDDVAR var = evbddnode_getvar(node);

    // Recursive apply to children and makenode (TODO: lace?)
    low  = evbdd_replace_terminal(low, b);
    high = evbdd_replace_terminal(high, b);
    res = evbdd_makenode(var, low, high);
    
    // We assume the input EVBDDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(EVBDD_WEIGHT(res) == EVBDD_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_EVBDD_REPLACE_TERMINAL, EVBDD_TARGET(a), b, 0, EVBDD_TARGET(res));
    return evbdd_bundle(EVBDD_TARGET(res), EVBDD_WEIGHT(a));
}

EVBDD
evbdd_tensor_prod(EVBDD a, EVBDD b, BDDVAR nvars_a)
{
    // Shift all vars of 'b' by 'nvars_a'
    b = evbdd_increase_all_vars(b, nvars_a);
 
    // Stack 'a' on top of 'b' (and multiply root edge of 'a' with that of 'b')
    EVBDD res = evbdd_replace_terminal(a, EVBDD_TARGET(b));
    EVBDD_WGT new_weight = wgt_mul(EVBDD_WEIGHT(res), EVBDD_WEIGHT(b));
    res = evbdd_bundle(EVBDD_TARGET(res), new_weight);

    return res;
}

/*************************</Matrix/vector operations>**************************/





/***************************<EVBDD utility functions>***************************/


/**
 * Helper function for recursive unmarking
 */
static void
evbdd_unmark_rec(EVBDD a)
{
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) return;
    evbddnode_t n = EVBDD_GETNODE(EVBDD_TARGET(a));
    if (!evbddnode_getmark(n)) return;
    evbddnode_setmark(n, 0);
    evbdd_unmark_rec(evbddnode_getptrlow(n));
    evbdd_unmark_rec(evbddnode_getptrhigh(n));
}

/**
 * Counts nodes in the EVBDD by marking them.
 */
static uint64_t
evbdd_nodecount_mark(EVBDD a)
{
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) return 0; // don't (repeat) count terminal
    evbddnode_t n = EVBDD_GETNODE(EVBDD_TARGET(a));
    if (evbddnode_getmark(n)) return 0;
    evbddnode_setmark(n, 1);
    return 1 + evbdd_nodecount_mark(evbddnode_getptrlow(n)) + evbdd_nodecount_mark(evbddnode_getptrhigh(n));
}

uint64_t
evbdd_countnodes(EVBDD a)
{
    uint64_t res = evbdd_nodecount_mark(a) + 1; // (+ 1 for terminal "node")
    evbdd_unmark_rec(a);
    return res;
}

/**************************</EVBDD utility functions>***************************/





/**************************<Printing & file writing>***************************/

/**
 * Pretty prints the information contained in `n`.
 */
static void evbddnode_pprint(evbddnode_t n)
{
    BDDVAR var = evbddnode_getvar(n);
    EVBDD low, high;
    evbddnode_getchilderen(n, &low, &high);
    printf("[var=%d, low=%" PRIu64 ", high=%" PRIu64 ", ",
             var,
             EVBDD_TARGET(low),
             EVBDD_TARGET(high));
    if(EVBDD_WEIGHT(low) == EVBDD_ZERO)      printf("a=EVBDD_ZERO, ");
    else if(EVBDD_WEIGHT(high) == EVBDD_ONE)  printf("a=EVBDD_ONE, ");
    else {
        printf("a=%" PRIu64 ", ", EVBDD_WEIGHT(low));
        printf("("); wgt_fprint(stdout, EVBDD_WEIGHT(low)); printf(")");
    }                      
    if(EVBDD_WEIGHT(high) == EVBDD_ZERO)     printf("b=EVBDD_ZERO ");
    else if(EVBDD_WEIGHT(high) == EVBDD_ONE) printf("b=EVBDD_ONE, ");
    else {                     
        printf("b=%" PRIu64, EVBDD_WEIGHT(high));
        printf("("); wgt_fprint(stdout, EVBDD_WEIGHT(high)); printf(")");
    }
    printf("]\n");
}

void
_print_evbdd(EVBDD a)
{
    if(EVBDD_TARGET(a) != EVBDD_TERMINAL){
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
        if(!evbddnode_getmark(node)){
            evbddnode_setmark(node, 1);
            printf("%" PRIu64 "\t", EVBDD_TARGET(a));
            evbddnode_pprint(node);
            EVBDD low, high;
            evbddnode_getchilderen(node, &low, &high);
            _print_evbdd(low);
            _print_evbdd(high);
        }
    }
}

void
evbdd_printnodes(EVBDD a)
{
    printf("root edge: %" PRIu64 ", %" PRIu64 " = ",EVBDD_TARGET(a), EVBDD_WEIGHT(a));
    wgt_fprint(stdout, EVBDD_WEIGHT(a));
    printf("\n");
    _print_evbdd(a);
    evbdd_unmark_rec(a);
}

static void
evbdd_fprintdot_edge_label(FILE *out, EVBDD_WGT w)
{
    fprintf(out, ", label=\"");
    if (w == EVBDD_ONE) {}
    else if (w == EVBDD_ZERO) { fprintf(out, "0"); }
    else if (w == EVBDD_MIN_ONE) { fprintf(out, "-1"); }
    else { wgt_fprint(out, w); }
    fprintf(out, "\"");
}

static void
evbdd_fprintdot_rec(FILE *out, EVBDD a, bool draw_zeros)
{
    // terminal node
    if(EVBDD_TARGET(a) == EVBDD_TERMINAL) return;

    evbddnode_t n = EVBDD_GETNODE(EVBDD_TARGET(a));
    if (evbddnode_getmark(n)) return;
    evbddnode_setmark(n, 1);

    // add this node
    fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n",
            EVBDD_TARGET(a), evbddnode_getvar(n));

    
    // children of this node
    EVBDD low, high;
    evbddnode_getchilderen(n, &low, &high);
    evbdd_fprintdot_rec(out, low, draw_zeros);
    evbdd_fprintdot_rec(out, high, draw_zeros);

    // add edge from this node to each child (unless weight 0)
    if (draw_zeros || EVBDD_WEIGHT(low) != EVBDD_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed",
                    EVBDD_TARGET(a), EVBDD_TARGET(low));
        evbdd_fprintdot_edge_label(out, EVBDD_WEIGHT(low));
        fprintf(out, "];\n");
    }
    if (draw_zeros || EVBDD_WEIGHT(high) != EVBDD_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                    EVBDD_TARGET(a), EVBDD_TARGET(high));
        evbdd_fprintdot_edge_label(out, EVBDD_WEIGHT(high));
        fprintf(out, "];\n");
    }
}

void
evbdd_fprintdot(FILE *out, EVBDD a, bool draw_zeros)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir=forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid", EVBDD_TARGET(a));
    evbdd_fprintdot_edge_label(out, EVBDD_WEIGHT(a));
    fprintf(out, "];\n");

    // terminal node
    fprintf(out, "%" PRIu64 " [shape=box, label=\"T\"];\n", EVBDD_TERMINAL);

    // recursively add nodes
    evbdd_fprintdot_rec(out, a, draw_zeros);
    evbdd_unmark_rec(a);

    fprintf(out, "}\n");
}

/**************************</Printing & file writing>**************************/





/********************************<Debug stuff>*********************************/

bool
evbdd_equivalent(EVBDD a, EVBDD b, int n, bool exact, bool verbose)
{
    bool has_next = true;
    EVBDD_WGT wgt_a, wgt_b;
    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    while(has_next){
        wgt_a = evbdd_getvalue(a, x);
        wgt_b = evbdd_getvalue(b, x);
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
evbdd_is_ordered_rec(EVBDD a, BDD nextvar, BDD nvars)
{
    // Terminal case
    if (EVBDD_TARGET(a) == EVBDD_TERMINAL) return true;

    // Get top node
    BDDVAR var;
    EVBDD low, high;
    evbdd_get_topvar(a, EVBDD_INVALID_VAR, &var, &low, &high);

    // Check variable order
    if (var >= nvars || var < nextvar) return false;

    // Check cache
    uint64_t res;
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_EVBDD_IS_ORDERED, EVBDD_TARGET(a), 0, 0, &res)) {
            return (bool) res;
        }
    }

    // Recursive calls
    bool res_low, res_high;
    var++;
    res_low = evbdd_is_ordered_rec(low, var, nvars);
    res_high = evbdd_is_ordered_rec(high, var, nvars);
    var--;
    res = (res_low && res_high);

    // Put res in cache and return
    if (cachenow) cache_put3(CACHE_EVBDD_IS_ORDERED, EVBDD_TARGET(a), 0, 0, res);
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
