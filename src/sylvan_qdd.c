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

#include "sylvan_int.h"
#include "sylvan_qdd.h"
#include "sylvan_refs.h"
#include "sylvan_qdd_complex.h"

static int granularity = 1; // operation cache access granularity
static bool testing_mode = 0; // turns on/off (expensive) sanity checks
static bool larger_amp_indices = false; // using [amps,ptr] [33,30] bits (default [23,40])
static int weight_norm_strat = 0;

static AMP (*normalize_weights)(AMP *, AMP *);

/****************< (bit level) manipulation of QDD / qddnode_t >***************/
/**
 * QDD edge structure (64 bits)
 *       1 bit:  marked/unmarked flag (same place as MTBDD)
 *      33 bits: index of edge weight in ctable (AMP)
 *      30 bits: index of next node in node table (PTR)
 * 
 * QDD node structure (128 bits)
 * (note: because of normalization of the amps, we only need 1 amp per node,
 *  the other will always be 0 or 1)
 * 
 * 64 bits low:
 *       1 bit:  marked/unmarked flag (same place as MTBDD)
 *       8 bits: variable/qubit number of this node
 *       1 bit:  if 0 (1) normalized amp is on low (high)
 *       1 bit:  if 0 (1) normalized amp is C_ZERO (C_ONE)
 *      13 bits: unused
 *      30 bits: low edge pointer to next node (PTR)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag (same place as MTBDD)
 *      33 bits: index of edge weight of high edge in ctable (AMP)
 *      30 bits: high edge pointer to next node (PTR)
 */
static const QDD qdd_marked_mask  = 0x8000000000000000LL;
static const QDD qdd_var_mask_low = 0x7f80000000000000LL;
static const QDD qdd_amp_pos_mask = 0x0040000000000000LL;
static const QDD qdd_amp_val_mask = 0x0020000000000000LL;
static const QDD qdd_amp_mask_23  = 0x7fffff0000000000LL;
static const QDD qdd_amp_mask_33  = 0x7fffffffc0000000LL;
static const QDD qdd_ptr_mask_30  = 0x000000003fffffffLL;
static const QDD qdd_ptr_mask_40  = 0x000000ffffffffffLL;

/**
 * Gets only the AMP information of a QDD edge `q`.
 */
static inline AMP
QDD_AMP(QDD q)
{
    if (larger_amp_indices) {
        return (q & qdd_amp_mask_33) >> 30; // 33 bits
    }
    else {
        return (q & qdd_amp_mask_23) >> 40; // 23 bits
    }
}

/**
 * Gets only the PTR information of a QDD edge `q`.
 */
static inline PTR
QDD_PTR(QDD q)
{
    if (larger_amp_indices) {
        return q & qdd_ptr_mask_30; // 30 bits
    }
    else {
        return q & qdd_ptr_mask_40; // 40 bits
    }
}

/**
 * 2 BDDVARs (assumed max 8 bits each)
 */
static inline uint32_t
QDD_PARAM_PACK_16(BDDVAR a, BDDVAR b) 
{
    return b<<8 | a;
}

/**
 *  24 bit gateid, 2 possible qubit parameters (e.g. control/target)
 */
static inline uint64_t
GATE_OPID_40(uint32_t gateid, BDDVAR a, BDDVAR b)
{
    uint64_t res = ((uint64_t)b)<<32 | ((uint64_t)a)<<24 | gateid;
    return res;
}

/**
 * 24 bits gateid, 5 possible qubit parameters (e.g. target/control/range)
 */
static inline uint64_t
GATE_OPID_64(uint32_t gateid, BDDVAR a, BDDVAR b, BDDVAR c, BDDVAR d, BDDVAR e)
{
    uint64_t res = ((uint64_t)e)<<56 | 
                   ((uint64_t)d)<<48 | 
                   ((uint64_t)c)<<40 | 
                   ((uint64_t)b)<<32 | 
                   ((uint64_t)a)<<24 | 
                   gateid;
    return res;
}

/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
qddnode_getvar(qddnode_t n)
{
    return (BDDVAR) ((n->low & qdd_var_mask_low) >> 55 ); // 8 bits
}

/**
 * Gets only the PTR of the low edge of `n`.
 */
static inline PTR
qddnode_getptrlow(qddnode_t n)
{
    return (PTR) QDD_PTR(n->low);
}

/**
 * Gets only the PTR of the high edge of `n`.
 */
static inline PTR
qddnode_getptrhigh(qddnode_t n)
{
    return (PTR) QDD_PTR(n->high);
}

/**
 * Gets the value of the "marked" flag.
 */
static inline bool
qddnode_getmark(qddnode_t n)
{
    return n->high & qdd_marked_mask ? 1 : 0;
}

/**
 * Sets the value of the "marked" flag to `mark`.
 */
static inline void
qddnode_setmark(qddnode_t n, bool mark)
{
    if (mark) n->high |=  qdd_marked_mask; // set 1st bit from left to 1
    else      n->high &= ~qdd_marked_mask; // set 1st bit from left to 0
}

/**
 * Gets the node `p` is pointing to.
 * TODO (?) return special node for when p == QDD_TERMINAL
 */
static inline qddnode_t
QDD_GETNODE(PTR p)
{
    return (qddnode_t) llmsset_index_to_ptr(nodes, p);
}

/**
 * Packs a PTR and AMP into a single 64 bit QDD.
 */
static inline QDD
qdd_bundle_ptr_amp(PTR p, AMP a)
{
    if (larger_amp_indices) {
        assert (p <= 0x000000003ffffffe);   // avoid clash with sylvan_invalid
        assert (a <= (1LL<<33));
        return (a << 30 | p);
    }else {
        assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
        assert (a <= (1<<23));
        return (a << 40 | p);
    }
}

static void
qddnode_unpack(qddnode_t n, PTR *low, PTR *high, AMP *a, AMP *b)
{
    *low  = qddnode_getptrlow(n);
    *high = qddnode_getptrhigh(n);
    bool norm_pos = (n->low & qdd_amp_pos_mask) >> 54;
    bool norm_val = (n->low & qdd_amp_val_mask) >> 53;

    if (weight_norm_strat == NORM_SUM) {
        *b = QDD_AMP(n->high);
        *a = amp_get_low_sum_normalized(*b);
    }
    else {
        if (norm_pos == 0) { // low amp is C_ZERO or C_ONE, high amp in ctable
            *a = (norm_val == 0) ? C_ZERO : C_ONE;
            *b = QDD_AMP(n->high);
        }
        else { // high amp is C_ZERO or C_ONE, low amp in ctable
            *b = (norm_val == 0) ? C_ZERO : C_ONE;
            *a = QDD_AMP(n->high);
        }
    }
}

static void
qddnode_getchilderen(qddnode_t n, QDD *low, QDD *high)
{
    PTR l, h;
    AMP a, b;
    qddnode_unpack(n, &l, &h, &a, &b);
    *low  = qdd_bundle_ptr_amp(l, a);
    *high = qdd_bundle_ptr_amp(h, b);
}

static void
qddnode_pack(qddnode_t n, BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    // We only want to store 1 complex number per node (which has 2 outgoing
    // edges). For NORM_LOW and NORM_LARGEST this is relatively easy because in
    // both those cases there is at least one edge weight equal to 1 or 0.
    //
    // For NORM_SUM it is a bit more complicated: both edge weights can be
    // outside of {0, 1}, but under the constraint that |low|^2 + |high|^2 = 1,
    // (or both are 0) and that |low| \in R+, we only need to store high, and
    // can derive low.

    // these will be set depending on the normalization strategy
    // (retrieval of edge weights is also dependent on normalization strategy)
    AMP amp_high;
    bool norm_pos;
    bool norm_val;
    
    if (weight_norm_strat == NORM_SUM) {
        assert(!(a == C_ZERO && b == C_ZERO)); // redundant node (caught before)
        norm_pos = 0;
        norm_val = 0;
        amp_high = b; // we can derive a from b
    }
    else {
        /// weight_norm_strat == NORM_LOW or NORM_LARGEST
        assert(a == C_ZERO || a == C_ONE || b == C_ZERO || b == C_ONE);
        norm_pos = (a == C_ZERO || a == C_ONE) ? 0 : 1;
        if (norm_pos == 0) {
            norm_val = (a == C_ZERO) ? 0 : 1;
            amp_high = b;
        }
        else {
            norm_val = (b == C_ZERO) ? 0 : 1;
            amp_high = a;
        }
    }

    // organize the bit structure of low and high
    n->low  = ((uint64_t)var)<<55 | ((uint64_t)norm_pos)<<54 | ((uint64_t)norm_val)<<53 | low;
    if (larger_amp_indices) {
        n->high = amp_high<<30 | high;
    }
    else {
        n->high = amp_high<<40 | high;
    }
}

// Container for disguising doubles as ints so they can go in Sylvan's cache
// (see also union "hack" in mtbdd_satcount)
typedef union {
    double   as_double;
    uint64_t as_int;
} double_hack_t;

/*
typedef union {
    complex_t as_comp;
    uint64_t  as_int[2];
} comp_hack_t;
*/

/***************</ (bit level) manipulation of QDD / qddnode_t >***************/



/*******************<garbage collection, references, marking>******************/

/**
 * Most of this gc code is copy-paste from sylvan_mtbdd.c, however because the 
 * bit structure of BDDs and QDDs are a bit different we can't use the mtbdd
 * code directly. Since the way sylvan_mtbdd/sylvan_common/sylvan_table are
 * structured we need to copy-paste a few more functions/variables than we 
 * actually change.
 */

/* 
 * Recursively mark QDD nodes as 'in use'.
 * This is really the only gc function which is different for QDDs vs MTBDDs.
 */
VOID_TASK_IMPL_1(qdd_gc_mark_rec, QDD, qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return;

    if (llmsset_mark(nodes, QDD_PTR(qdd))) {
        qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
        SPAWN(qdd_gc_mark_rec, qddnode_getptrlow(n));
        CALL(qdd_gc_mark_rec, qddnode_getptrhigh(n));
        SYNC(qdd_gc_mark_rec);
    }
}

/**
 * External references
 */
refs_table_t qdd_refs;
refs_table_t qdd_protected;
static int qdd_protected_created = 0;

void
qdd_protect(QDD *a)
{
    if (!qdd_protected_created) {
        // In C++, sometimes mtbdd_protect is called before Sylvan is initialized. Just create a table.
        protect_create(&qdd_protected, 4096);
        qdd_protected_created = 1;
    }
    protect_up(&qdd_protected, (size_t)a);
}

void
qdd_unprotect(QDD *a)
{
    if (qdd_protected.refs_table != NULL) protect_down(&qdd_protected, (size_t)a);
}

size_t
qdd_count_protected()
{
    return protect_count(&qdd_protected);
}

/* Called during garbage collection */
VOID_TASK_0(qdd_gc_mark_external_refs)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = refs_iter(&qdd_refs, 0, qdd_refs.refs_size);
    while (it != NULL) {
        SPAWN(qdd_gc_mark_rec, refs_next(&qdd_refs, &it, qdd_refs.refs_size));
        count++;
    }
    while (count--) {
        SYNC(qdd_gc_mark_rec);
    }
}

/* Called during garbage collection */
VOID_TASK_0(qdd_gc_mark_protected)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = protect_iter(&qdd_protected, 0, qdd_protected.refs_size);
    while (it != NULL) {
        QDD *to_mark = (QDD*)protect_next(&qdd_protected, &it, qdd_protected.refs_size);
        SPAWN(qdd_gc_mark_rec, *to_mark);
        count++;
    }
    while (count--) {
        SYNC(qdd_gc_mark_rec);
    }
}

/* Infrastructure for internal markings */
typedef struct qdd_refs_task
{
    Task *t;
    void *f;
} *qdd_refs_task_t;

typedef struct qdd_refs_internal
{
    const QDD **pbegin, **pend, **pcur;
    QDD *rbegin, *rend, *rcur;
    qdd_refs_task_t sbegin, send, scur;
} *qdd_refs_internal_t;

DECLARE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);

VOID_TASK_2(qdd_refs_mark_p_par, const QDD**, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            qdd_gc_mark_rec(**(begin++));
            count--;
        }
    } else {
        SPAWN(qdd_refs_mark_p_par, begin, count / 2);
        CALL(qdd_refs_mark_p_par, begin + (count / 2), count - count / 2);
        SYNC(qdd_refs_mark_p_par);
    }
}

VOID_TASK_2(qdd_refs_mark_r_par, QDD*, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            qdd_gc_mark_rec(*begin++);
            count--;
        }
    } else {
        SPAWN(qdd_refs_mark_r_par, begin, count / 2);
        CALL(qdd_refs_mark_r_par, begin + (count / 2), count - count / 2);
        SYNC(qdd_refs_mark_r_par);
    }
}

VOID_TASK_2(qdd_refs_mark_s_par, qdd_refs_task_t, begin, size_t, count)
{
    if (count < 32) {
        while (count > 0) {
            Task *t = begin->t;
            if (!TASK_IS_STOLEN(t)) return;
            if (t->f == begin->f && TASK_IS_COMPLETED(t)) {
                qdd_gc_mark_rec(*(QDD*)TASK_RESULT(t));
            }
            begin += 1;
            count -= 1;
        }
    } else {
        if (!TASK_IS_STOLEN(begin->t)) return;
        SPAWN(qdd_refs_mark_s_par, begin, count / 2);
        CALL(qdd_refs_mark_s_par, begin + (count / 2), count - count / 2);
        SYNC(qdd_refs_mark_s_par);
    }
}

VOID_TASK_0(qdd_refs_mark_task)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    SPAWN(qdd_refs_mark_p_par, qdd_refs_key->pbegin, qdd_refs_key->pcur-qdd_refs_key->pbegin);
    SPAWN(qdd_refs_mark_r_par, qdd_refs_key->rbegin, qdd_refs_key->rcur-qdd_refs_key->rbegin);
    CALL(qdd_refs_mark_s_par, qdd_refs_key->sbegin, qdd_refs_key->scur-qdd_refs_key->sbegin);
    SYNC(qdd_refs_mark_r_par);
    SYNC(qdd_refs_mark_p_par);
}

/* Called during garbage collection */
VOID_TASK_0(qdd_refs_mark)
{
    TOGETHER(qdd_refs_mark_task);
}

VOID_TASK_0(qdd_refs_init_task)
{
    qdd_refs_internal_t s = (qdd_refs_internal_t)malloc(sizeof(struct qdd_refs_internal));
    s->pcur = s->pbegin = (const QDD**)malloc(sizeof(QDD*) * 1024);
    s->pend = s->pbegin + 1024;
    s->rcur = s->rbegin = (QDD*)malloc(sizeof(QDD) * 1024);
    s->rend = s->rbegin + 1024;
    s->scur = s->sbegin = (qdd_refs_task_t)malloc(sizeof(struct qdd_refs_task) * 1024);
    s->send = s->sbegin + 1024;
    SET_THREAD_LOCAL(qdd_refs_key, s);
}

VOID_TASK_0(qdd_refs_init)
{
    INIT_THREAD_LOCAL(qdd_refs_key);
    TOGETHER(qdd_refs_init_task);
    sylvan_gc_add_mark(TASK(qdd_refs_mark));
}

/**
 * Initialize and quit functions
 */
static int qdd_initialized = 0;

static void
qdd_quit()
{
    refs_free(&qdd_refs);
    if (qdd_protected_created) {
        protect_free(&qdd_protected);
        qdd_protected_created = 0;
    }
    qdd_initialized = 0;
    free_amplitude_table();
}

void
sylvan_init_qdd(size_t ctable_size, double ctable_tolerance, int amps_backend, int norm_strat)
{
    if (qdd_initialized) return;
    qdd_initialized = 1;

    int index_size = (int) ceil(log2(ctable_size));
    if (index_size > 33) {
        printf("max amp storage size is 2^33 (half when storing real components separately)\n");
        exit(1);
    }
    if (index_size > 23) {
        larger_amp_indices = true;
    }

    sylvan_register_quit(qdd_quit);
    sylvan_gc_add_mark(TASK(qdd_gc_mark_external_refs));
    sylvan_gc_add_mark(TASK(qdd_gc_mark_protected));

    refs_create(&qdd_refs, 1024);
    if (!qdd_protected_created) {
        protect_create(&qdd_protected, 4096);
        qdd_protected_created = 1;
    }

    init_amplitude_table(ctable_size, ctable_tolerance, amps_backend);
    qdd_gates_init();

    weight_norm_strat = norm_strat;
    switch (norm_strat)
    {
    case NORM_LARGEST:
        normalize_weights = amp_normalize_largest_ptr;
        break;
    case NORM_LOW:
        normalize_weights = amp_normalize_low_ptr;
        break;
    case NORM_SUM:
        normalize_weights = amp_normalize_sum_ptr;
        break;
    default:
        printf("Edge weight normalization strategy not recognized\n");
        exit(1);
        break;
    }

    LACE_ME;
    CALL(qdd_refs_init);
}

void
sylvan_init_qdd_defaults(size_t ctable_size)
{
    sylvan_init_qdd(ctable_size, -1, COMP_HASHMAP, NORM_LOW);
}

void
qdd_set_testing_mode(bool on)
{
    testing_mode = on;
}

void
qdd_set_caching_granularity(int g)
{
    granularity = g;
}

void
qdd_refs_ptrs_up(qdd_refs_internal_t qdd_refs_key)
{
    size_t cur = qdd_refs_key->pcur - qdd_refs_key->pbegin;
    size_t size = qdd_refs_key->pend - qdd_refs_key->pbegin;
    qdd_refs_key->pbegin = (const QDD**)realloc(qdd_refs_key->pbegin, sizeof(QDD*) * size * 2);
    qdd_refs_key->pcur = qdd_refs_key->pbegin + cur;
    qdd_refs_key->pend = qdd_refs_key->pbegin + (size * 2);
}

QDD __attribute__((noinline))
qdd_refs_refs_up(qdd_refs_internal_t qdd_refs_key, QDD res)
{
    long size = qdd_refs_key->rend - qdd_refs_key->rbegin;
    qdd_refs_key->rbegin = (QDD*)realloc(qdd_refs_key->rbegin, sizeof(QDD) * size * 2);
    qdd_refs_key->rcur = qdd_refs_key->rbegin + size;
    qdd_refs_key->rend = qdd_refs_key->rbegin + (size * 2);
    return res;
}

void __attribute__((noinline))
qdd_refs_tasks_up(qdd_refs_internal_t qdd_refs_key)
{
    long size = qdd_refs_key->send - qdd_refs_key->sbegin;
    qdd_refs_key->sbegin = (qdd_refs_task_t)realloc(qdd_refs_key->sbegin, sizeof(struct qdd_refs_task) * size * 2);
    qdd_refs_key->scur = qdd_refs_key->sbegin + size;
    qdd_refs_key->send = qdd_refs_key->sbegin + (size * 2);
}

void __attribute__((unused))
qdd_refs_pushptr(const QDD *ptr)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    *qdd_refs_key->pcur++ = ptr;
    if (qdd_refs_key->pcur == qdd_refs_key->pend) qdd_refs_ptrs_up(qdd_refs_key);
}

void __attribute__((unused))
qdd_refs_popptr(size_t amount)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    qdd_refs_key->pcur -= amount;
}

QDD __attribute__((unused))
qdd_refs_push(QDD qdd)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    *(qdd_refs_key->rcur++) = qdd;
    if (qdd_refs_key->rcur == qdd_refs_key->rend) return qdd_refs_refs_up(qdd_refs_key, qdd);
    else return qdd;
}

void __attribute__((unused))
qdd_refs_pop(long amount)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    qdd_refs_key->rcur -= amount;
}

void
qdd_refs_spawn(Task *t)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    qdd_refs_key->scur->t = t;
    qdd_refs_key->scur->f = t->f;
    qdd_refs_key->scur += 1;
    if (qdd_refs_key->scur == qdd_refs_key->send) qdd_refs_tasks_up(qdd_refs_key);
}

QDD
qdd_refs_sync(QDD result)
{
    LOCALIZE_THREAD_LOCAL(qdd_refs_key, qdd_refs_internal_t);
    qdd_refs_key->scur -= 1;
    return result;
}

/******************</garbage collection, references, marking>******************/

/**
 * Gets either the top node of the qdd, or the node with var 't' if this
 * variable would otherwise be skipped. The "node" is returned as 
 * (*topvar, *low, *high).
 */
static void
qdd_get_topvar(QDD qdd, BDDVAR t, BDDVAR *topvar, QDD *low, QDD *high)
{
    bool skipped = false;
    if(QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
        *topvar = qddnode_getvar(node);
        if (*topvar > t) skipped = true;
    }

    if (skipped) {
        *low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        *high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        *topvar = t;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
        qddnode_getchilderen(node, low, high);
    }
}

static PTR
_qdd_makenode(BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    struct qddnode n;

    /*
    int mark;
    if (MTBDD_HASMARK(low)) {
        mark = 1;
        low = MTBDD_TOGGLEMARK(low);
        high = MTBDD_TOGGLEMARK(high);
    } else {
        mark = 0;
    }
    */

    qddnode_pack(&n, var, low, high, a, b);

    PTR result;
    int created;
    PTR index = llmsset_lookup(nodes, n.low, n.high, &created);
    if (index == 0) {
        printf("auto gc of node table triggered\n");
        LACE_ME;

        qdd_refs_push(low);
        qdd_refs_push(high);
        sylvan_gc();
        qdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.low, n.high, &created);
        if (index == 0) {
            fprintf(stderr, "QDD/BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(QDD_NODES_CREATED);
    else sylvan_stats_count(QDD_NODES_REUSED);

    result = index;
    //return mark ? result | qdd_marked_mask : result;
    return result;
}


static AMP
qdd_comp_lookup(complex_t c)
{
    // !! TODO: replace with try_lookup, catch failure, trigger gc amp table
    AMP res;
    bool inserted;
    res = comp_try_lookup(c, &inserted);
    if (!inserted) {
        printf("TODO: GC AMP TABLE\n");
        exit(1);
    }
    return res;
}

static QDD // (PTR and AMP, but the amp is the norm weight from below)
qdd_makenode(BDDVAR var, QDD low, QDD high)
{ 
    PTR low_ptr  = QDD_PTR(low);
    AMP low_amp  = QDD_AMP(low);
    PTR high_ptr = QDD_PTR(high);
    AMP high_amp = QDD_AMP(high);

    // Edges with weight 0 should point straight to terminal.
    if (low_amp  == C_ZERO) low_ptr  = QDD_TERMINAL;
    if (high_amp == C_ZERO) high_ptr = QDD_TERMINAL;

    // If both low and high are the same (both PTR and AMP) return low
    if (low == high) return low;
    else {
        // If the edges are not the same
        AMP norm = (*normalize_weights)(&low_amp, &high_amp);
        PTR res  = _qdd_makenode(var, low_ptr, high_ptr, low_amp, high_amp);
        return qdd_bundle_ptr_amp(res, norm);
    }
}


/**
 * Helper function for recursive unmarking
 */
static void
qdd_unmark_rec(QDD qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return;
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (!qddnode_getmark(n)) return;
    qddnode_setmark(n, 0);
    qdd_unmark_rec(qddnode_getptrlow(n));
    qdd_unmark_rec(qddnode_getptrhigh(n));
}

/**
 * Counts nodes in the qdd by marking them.
 */
static uint64_t
qdd_nodecount_mark(QDD qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return 0; // don't (repeat) count terminal
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (qddnode_getmark(n)) return 0;
    qddnode_setmark(n, 1);
    return 1 + qdd_nodecount_mark(qddnode_getptrlow(n)) + qdd_nodecount_mark(qddnode_getptrhigh(n));
}

uint64_t
qdd_countnodes(QDD qdd)
{
    uint64_t res = qdd_nodecount_mark(qdd) + 1; // (+ 1 for terminal "node")
    qdd_unmark_rec(qdd);
    return res;
}


/**************************<cleaning amplitude table>**************************/

static int auto_gc_amp_table  = 1;
static double amp_table_gc_thres = 0.5;

void
qdd_set_auto_gc_amp_table(bool enabled)
{
    auto_gc_amp_table = enabled;
}

void
qdd_set_gc_amp_table_thres(double fraction_filled)
{
    amp_table_gc_thres = fraction_filled;
}

double
qdd_get_gc_amp_table_thres()
{
    return amp_table_gc_thres;
}


void
qdd_gc_amp_table()
{
    // gc amp table and keep amps of protected qdds (and update those)
    LACE_ME; 
    
    // 1. Create new amp table
    init_new_empty_table();

    // 2. Fill new table with amps in protected qdd's and update those qdd's
    uint64_t *it = protect_iter(&qdd_protected, 0, qdd_protected.refs_size);
    while (it != NULL) {
        QDD *to_protect_amps = (QDD*)protect_next(&qdd_protected, &it, qdd_protected.refs_size);
        if (to_protect_amps != NULL) {
            *to_protect_amps = _fill_new_amp_table(*to_protect_amps);
        }
    }

    // 3. Delete old amp table
    delete_old_table();

    // 4. Any cache we migh have is now invalid because the same amplitudes 
    //    might now have different indices in the amp table
    sylvan_clear_cache();
}

TASK_IMPL_1(QDD, _fill_new_amp_table, QDD, qdd)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CLEAN_AMP_TABLE, 0LL, qdd, 0LL, &res)) {
            return res;
        }
    }

    // Move amp from old to new table, get new index
    AMP new_amp = move_from_old_to_new(QDD_AMP(qdd));
    qdd = qdd_bundle_ptr_amp(QDD_PTR(qdd), new_amp);

    // If terminal, return
    if (QDD_PTR(qdd) == QDD_TERMINAL) return qdd;
    
    // Recursive for children
    QDD low, high;
    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    qddnode_getchilderen(n, &low, &high);
    qdd_refs_spawn(SPAWN(_fill_new_amp_table, high));
    low = CALL(_fill_new_amp_table, low);
    qdd_refs_push(low);
    high = qdd_refs_sync(SYNC(_fill_new_amp_table));
    qdd_refs_pop(1);

    // We don't need to use the 'qdd_makenode()' function which normalizes the 
    // amplitudes, because the QDD doesn't actually change, only the AMP indices,
    // but none of the actual values.
    PTR ptr = _qdd_makenode(qddnode_getvar(n), QDD_PTR(low), QDD_PTR(high), 
                                               QDD_AMP(low), QDD_AMP(high));

    // Put in cache, return
    res = qdd_bundle_ptr_amp(ptr, new_amp);
    if (cachenow) cache_put3(CACHE_QDD_CLEAN_AMP_TABLE, 0LL, qdd, 0LL, res);
    return res;
}

bool
qdd_test_gc_amp_table()
{
    uint64_t entries = get_table_entries_estimate();
    uint64_t size    = get_table_size();
    return ( ((double)entries / (double)size) > amp_table_gc_thres );
}

/*************************</cleaning amplitude table>**************************/



/*******************************<applying gates>*******************************/

// temp trigger for gc of node table
static int periodic_gc_nodetable = 0;
static uint64_t gate_counter = 0;
void
qdd_set_periodic_gc_nodetable(int every_n_gates)
{
    periodic_gc_nodetable = every_n_gates;
}

static void
qdd_do_before_gate(QDD* qdd)
{
    // check if ctable needs gc
    if (auto_gc_amp_table && qdd_test_gc_amp_table()) {
        qdd_protect(qdd);
        qdd_gc_amp_table();
        qdd_unprotect(qdd);
    }

    if (periodic_gc_nodetable) {
        gate_counter++;
        if (gate_counter % periodic_gc_nodetable ==  0) {
            LACE_ME;
            qdd_protect(qdd);
            sylvan_gc();
            qdd_unprotect(qdd);
        }
    }

    // log stuff (if logging is enabled)
    qdd_stats_log(*qdd);
}

static void
qdd_do_before_mult()
{
    // check if ctable needs gc
    if (auto_gc_amp_table && qdd_test_gc_amp_table()) {
        qdd_gc_amp_table();
    }
}

/* Wrapper for applying a single qubit gate. */
TASK_IMPL_3(QDD, qdd_gate, QDD, qdd, gate_id_t, gate, BDDVAR, target)
{
    qdd_do_before_gate(&qdd);
    return qdd_gate_rec(qdd, gate, target);
}

/* Wrapper for applying a controlled gate with 1 control qubit. */
TASK_IMPL_4(QDD, qdd_cgate, QDD, qdd, gate_id_t, gate, BDDVAR, c, BDDVAR, t)
{
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c, QDD_INVALID_VAR, QDD_INVALID_VAR, QDD_INVALID_VAR};
    return qdd_cgate_rec(qdd, gate, cs, t);
}

/* Wrapper for applying a controlled gate with 2 control qubits. */
TASK_IMPL_5(QDD, qdd_cgate2, QDD, qdd, gate_id_t, gate, BDDVAR, c1, BDDVAR, c2, BDDVAR, t)
{
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c1, c2, QDD_INVALID_VAR, QDD_INVALID_VAR};
    return qdd_cgate_rec(qdd, gate, cs, t);
}

/* Wrapper for applying a controlled gate with 3 control qubits. */
TASK_IMPL_6(QDD, qdd_cgate3, QDD, qdd, gate_id_t, gate, BDDVAR, c1, BDDVAR, c2, BDDVAR, c3, BDDVAR, t)
{
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c1, c2, c3, QDD_INVALID_VAR}; // last pos is a buffer
    return qdd_cgate_rec(qdd, gate, cs, t);
}

/* Wrapper for applying a controlled gate where the controls are a range. */
TASK_IMPL_5(QDD, qdd_cgate_range, QDD, qdd, gate_id_t, gate, BDDVAR, c_first, BDDVAR, c_last, BDDVAR, t)
{
    qdd_do_before_gate(&qdd);
    return qdd_cgate_range_rec(qdd,gate,c_first,c_last,t);
}

static void
norm_plus_cache_key(QDD a, QDD b, QDD *x, QDD *y)
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

TASK_IMPL_2(QDD, qdd_plus, QDD, a, QDD, b)
{
    // Trivial cases
    if(QDD_AMP(a) == C_ZERO) return b;
    if(QDD_AMP(b) == C_ZERO) return a;

    // Get var(a) and var(b)
    QDD low_a, low_b, high_a, high_b, res;
    BDDVAR var_a = UINT32_MAX, var_b = UINT32_MAX, topvar;
    if (QDD_PTR(a) != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(QDD_PTR(a));
        var_a  = qddnode_getvar(node);
    }
    if (QDD_PTR(b) != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(QDD_PTR(b));
        var_b  = qddnode_getvar(node);
    }

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    qdd_get_topvar(a, var_b, &topvar, &low_a, &high_a);
    qdd_get_topvar(b, var_a, &topvar, &low_b, &high_b);

    // Base/terminal case: same target and same variable
    if(QDD_PTR(a) == QDD_PTR(b) && var_a == var_b){
        AMP sum = amp_add(QDD_AMP(a), QDD_AMP(b));
        res = qdd_bundle_ptr_amp(QDD_PTR(a), sum);
        return res;
    }

    // Check cache
    QDD x, y;
    norm_plus_cache_key(a, b, &x, &y); // (a + b) = (b + a) so normalize cache key
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_PLUS, sylvan_false, x, y, &res)) {
            sylvan_stats_count(QDD_PLUS_CACHED);
            return res;
        }
    }

    // If not base/terminal case, pass edge weight of current edge down
    AMP amp_la, amp_ha, amp_lb, amp_hb;
    amp_la = amp_mul_down(QDD_AMP(a), QDD_AMP(low_a));
    amp_ha = amp_mul_down(QDD_AMP(a), QDD_AMP(high_a));
    amp_lb = amp_mul_down(QDD_AMP(b), QDD_AMP(low_b));
    amp_hb = amp_mul_down(QDD_AMP(b), QDD_AMP(high_b));
    low_a  = qdd_bundle_ptr_amp(QDD_PTR(low_a),  amp_la);
    high_a = qdd_bundle_ptr_amp(QDD_PTR(high_a), amp_ha);
    low_b  = qdd_bundle_ptr_amp(QDD_PTR(low_b),  amp_lb);
    high_b = qdd_bundle_ptr_amp(QDD_PTR(high_b), amp_hb);

    // Recursive calls down
    QDD low, high;
    qdd_refs_spawn(SPAWN(qdd_plus, high_a, high_b));
    low = CALL(qdd_plus, low_a, low_b);
    qdd_refs_push(low);
    high = qdd_refs_sync(SYNC(qdd_plus));
    qdd_refs_pop(1);

    // Put in cache, return
    res = qdd_makenode(topvar, low, high);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_PLUS, sylvan_false, x, y, res)) 
            sylvan_stats_count(QDD_PLUS_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_3(QDD, qdd_gate_rec, QDD, q, gate_id_t, gate, BDDVAR, target)
{
    // Trivial cases
    if (QDD_AMP(q) == C_ZERO) return q;

    BDDVAR var;
    QDD res, low, high;
    qdd_get_topvar(q, target, &var, &low, &high);
    assert(var <= target);

    // Check cache
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_GATE, sylvan_false, QDD_PTR(q), GATE_OPID_40(gate, target, 0), &res)) {
            sylvan_stats_count(QDD_GATE_CACHED);
            // Multiply root of res with root of input qdd
            AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
            res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
            return res;
        }
    }

    if (var == target) {
        AMP a_u00 = amp_mul(QDD_AMP(low), gates[gate][0]);
        AMP a_u10 = amp_mul(QDD_AMP(low), gates[gate][2]);
        AMP b_u01 = amp_mul(QDD_AMP(high), gates[gate][1]);
        AMP b_u11 = amp_mul(QDD_AMP(high), gates[gate][3]);
        QDD low1, low2, high1, high2;
        low1  = qdd_bundle_ptr_amp(QDD_PTR(low), a_u00);
        low2  = qdd_bundle_ptr_amp(QDD_PTR(high),b_u01);
        high1 = qdd_bundle_ptr_amp(QDD_PTR(low), a_u10);
        high2 = qdd_bundle_ptr_amp(QDD_PTR(high),b_u11);
        qdd_refs_spawn(SPAWN(qdd_plus, high1, high2));
        low = CALL(qdd_plus, low1, low2);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_plus));
        qdd_refs_pop(1);
        res = qdd_makenode(target, low, high);
    }
    else { // var < target: not at target qubit yet, recursive calls down
        qdd_refs_spawn(SPAWN(qdd_gate_rec, high, gate, target));
        low = CALL(qdd_gate_rec, low, gate, target);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_gate_rec));
        qdd_refs_pop(1);
        res  = qdd_makenode(var, low, high);
    }

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_GATE, sylvan_false, QDD_PTR(q), GATE_OPID_40(gate, target, 0), res)) 
            sylvan_stats_count(QDD_GATE_CACHEDPUT);
    }
    // Multiply amp res with amp of input qdd
    AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

TASK_IMPL_5(QDD, qdd_cgate_rec, QDD, q, gate_id_t, gate, BDDVAR*, cs, uint32_t, ci, BDDVAR, t)
{
    // Get current control qubit. If no more control qubits, apply gate here
    BDDVAR c = cs[ci];
    if (c == QDD_INVALID_VAR || ci > MAX_CONTROLS) {
        return CALL(qdd_gate_rec, q, gate, t);
    }

    assert(c < t && "ctrl < target required");
    if (ci > 0) 
        assert(cs[ci-1] < cs[ci]  && "order required for multiple controls");

    BDDVAR var;
    QDD res, low, high;
    qdd_get_topvar(q, c, &var, &low, &high);
    assert(var <= c);

    // Check cache
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CGATE, sylvan_false, QDD_PTR(q),
                    GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                    &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            // Multiply root amp of res with input root amp
            AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
            res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
            return res;
        }
    }

    // If current node is (one of) the control qubit(s), 
    // control on q_c = |1> (high edge)
    if (var == c) {
        ci++;
        high = CALL(qdd_cgate_rec, high, gate, cs, ci, t);
        ci--;
    }
    // Not at control qubit yet, apply to both childeren.
    else {
        qdd_refs_spawn(SPAWN(qdd_cgate_rec, high, gate, cs, ci, t));
        low = CALL(qdd_cgate_rec, low, gate, cs, ci, t);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_cgate_rec));
        qdd_refs_pop(1);
    }
    res = qdd_makenode(var, low, high);

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE, sylvan_false, QDD_PTR(q),
                    GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                    res)) {
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
        }
    }
    // Multiply root amp of res with input root amp
    AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

TASK_IMPL_6(QDD, qdd_cgate_range_rec, QDD, q, gate_id_t, gate, BDDVAR, c_first, BDDVAR, c_last, BDDVAR, t, BDDVAR, k)
{
    // Past last control (done with "control part" of controlled gate)
    if (k > c_last) {
        return CALL(qdd_gate_rec, q, gate, t);
    }

    assert(c_first <= c_last);
    assert(c_last < t);

    // Check cache
    QDD res;
    bool cachenow = ((k % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CGATE_RANGE, sylvan_false, QDD_PTR(q),
                       GATE_OPID_64(gate, c_first, c_last, k, t, 0),
                       &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            // Multiply root amp of result with the input root amp
            AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
            res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
            return res;
        }
    }

    // Get top node
    BDDVAR var, nextvar;
    QDD low, high;
    if (k < c_first) {
        // possibly skip to c_first if we are before c_first
        qdd_get_topvar(q, c_first, &var, &low, &high);
        assert(var <= c_first);
    }
    else {
        // k is a control qubit, so get node with var = k (insert if skipped)
        assert(c_first <= k && k <= c_last);
        qdd_get_topvar(q, k, &var, &low, &high);
        assert(var == k);
    }
    nextvar = var + 1;
    
    // Not at first control qubit yet, apply to both children
    if (var < c_first) {
        qdd_refs_spawn(SPAWN(qdd_cgate_range_rec, high, gate, c_first, c_last, t, nextvar));
        low = CALL(qdd_cgate_range_rec, low, gate, c_first, c_last, t, var);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_cgate_range_rec));
        qdd_refs_pop(1);
    }
    // Current var is a control qubit, control on q_k = |1> (high edge)
    else {
        high = CALL(qdd_cgate_range_rec, high, gate, c_first, c_last, t, nextvar);
    }
    res = qdd_makenode(var, low, high);

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE_RANGE, sylvan_false, QDD_PTR(q),
                       GATE_OPID_64(gate, c_first, c_last, k, t, 0),
                       res)) {
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
        }
    }
    // Multiply root amp of result with the input root amp
    AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(QDD, qdd_matvec_mult, QDD, mat, QDD, vec, BDDVAR, nvars)
{
    qdd_do_before_mult();
    return CALL(qdd_matvec_mult_rec, mat, vec, nvars, 0);
}

/* Wrapper for matrix vector multiplication. */
TASK_IMPL_3(QDD, qdd_matmat_mult, QDD, a, QDD, b, BDDVAR, nvars)
{
    qdd_do_before_mult();
    return CALL(qdd_matmat_mult_rec, a, b, nvars, 0);
}

TASK_IMPL_4(QDD, qdd_matvec_mult_rec, QDD, mat, QDD, vec, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (QDD_AMP(mat) == C_ZERO || QDD_AMP(vec) == C_ZERO)
        return qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
    
    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(QDD_PTR(mat) == QDD_TERMINAL);
        assert(QDD_PTR(vec) == QDD_TERMINAL);
        AMP prod = amp_mul(QDD_AMP(mat), QDD_AMP(vec));
        return qdd_bundle_ptr_amp(QDD_TERMINAL, prod);
    }

    // Check cache
    QDD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_MATVEC_MULT, nextvar, QDD_PTR(mat), QDD_PTR(vec), &res)) {
            sylvan_stats_count(QDD_MULT_CACHED);
            // 6. multiply w/ product of root amps
            AMP prod = amp_mul(QDD_AMP(mat), QDD_AMP(vec));
            AMP new_root_amp = amp_mul(prod, QDD_AMP(res));
            res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both QDDS
    BDDVAR var;
    QDD vec_low, vec_high, mat_low, mat_high, u00, u10, u01, u11;
    qdd_get_topvar(vec, nextvar, &var, &vec_low, &vec_high);
    qdd_get_topvar(mat, 2*nextvar, &var, &mat_low, &mat_high);
    qdd_get_topvar(mat_low, 2*nextvar+1, &var, &u00, &u10);
    qdd_get_topvar(mat_high,2*nextvar+1, &var, &u01, &u11);

    // 2. propagate "in-between" amps of matrix qdd
    u00 = qdd_bundle_ptr_amp(QDD_PTR(u00), amp_mul(QDD_AMP(u00), QDD_AMP(mat_low)));
    u10 = qdd_bundle_ptr_amp(QDD_PTR(u10), amp_mul(QDD_AMP(u10), QDD_AMP(mat_low)));
    u01 = qdd_bundle_ptr_amp(QDD_PTR(u01), amp_mul(QDD_AMP(u01), QDD_AMP(mat_high)));
    u11 = qdd_bundle_ptr_amp(QDD_PTR(u11), amp_mul(QDD_AMP(u11), QDD_AMP(mat_high)));

    // 3. recursive calls (4 tasks: SPAWN 3, CALL 1)
    // |u00 u01| |vec_low | = vec_low|u00| + vec_high|u01|
    // |u10 u11| |vec_high|          |u10|           |u11|
    QDD res_low00, res_low10, res_high01, res_high11;
    nextvar++;
    qdd_refs_spawn(SPAWN(qdd_matvec_mult_rec, u00, vec_low,  nvars, nextvar)); // 1
    qdd_refs_spawn(SPAWN(qdd_matvec_mult_rec, u10, vec_low,  nvars, nextvar)); // 2
    qdd_refs_spawn(SPAWN(qdd_matvec_mult_rec, u01, vec_high, nvars, nextvar)); // 3
    res_high11 = CALL(qdd_matvec_mult_rec, u11, vec_high, nvars, nextvar);
    qdd_refs_push(res_high11);
    res_high01 = qdd_refs_sync(SYNC(qdd_matvec_mult_rec)); // 3
    res_low10  = qdd_refs_sync(SYNC(qdd_matvec_mult_rec)); // 2
    res_low00  = qdd_refs_sync(SYNC(qdd_matvec_mult_rec)); // 1
    qdd_refs_pop(1);
    nextvar--;

    // 4. gather results of multiplication
    QDD res_low, res_high;
    res_low  = qdd_makenode(nextvar, res_low00,  res_low10);
    res_high = qdd_makenode(nextvar, res_high01, res_high11);

    // 5. add resulting qdds
    res = CALL(qdd_plus, res_low, res_high);

    // Insert in cache (before multiplication w/ root amps)
    if (cachenow) {
        if (cache_put3(CACHE_QDD_MATVEC_MULT, nextvar, QDD_PTR(mat), QDD_PTR(vec), res)) 
            sylvan_stats_count(QDD_MULT_CACHEDPUT);
    }

    // 6. multiply w/ product of root amps
    AMP prod = amp_mul(QDD_AMP(mat), QDD_AMP(vec));
    AMP new_root_amp = amp_mul(prod, QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

    return res;
}

TASK_IMPL_4(QDD, qdd_matmat_mult_rec, QDD, a, QDD, b, BDDVAR, nvars, BDDVAR, nextvar)
{
    // Trivial case: either one is all 0
    if (QDD_AMP(a) == C_ZERO || QDD_AMP(b) == C_ZERO)
        return qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);

    // Terminal case: past last variable
    if (nextvar == nvars) {
        assert(QDD_PTR(a) == QDD_TERMINAL);
        assert(QDD_PTR(b) == QDD_TERMINAL);
        AMP prod = amp_mul(QDD_AMP(a), QDD_AMP(b));
        return qdd_bundle_ptr_amp(QDD_TERMINAL, prod);
    }

    // Check cache
    QDD res;
    bool cachenow = ((nextvar % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_MATMAT_MULT, nextvar, QDD_PTR(a), QDD_PTR(b), &res)) {
            sylvan_stats_count(QDD_MULT_CACHED);
            // 7. multiply w/ product of root amps
            AMP prod = amp_mul(QDD_AMP(a), QDD_AMP(b));
            AMP new_root_amp = amp_mul(prod, QDD_AMP(res));
            res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
            return res;
        }
    }

    // Recursive multiplication
    // 1. get relevant nodes for both QDDs
    BDDVAR var;
    QDD a_low, a_high, a00, a10, a01, a11, b_low, b_high, b00, b10, b01, b11;
    qdd_get_topvar(a, 2*nextvar, &var, &a_low, &a_high);
    qdd_get_topvar(b, 2*nextvar, &var, &b_low, &b_high);
    qdd_get_topvar(a_low, 2*nextvar+1, &var, &a00, &a10);
    qdd_get_topvar(a_high,2*nextvar+1, &var, &a01, &a11);
    qdd_get_topvar(b_low, 2*nextvar+1, &var, &b00, &b10);
    qdd_get_topvar(b_high,2*nextvar+1, &var, &b01, &b11);

    // 2. propagate "in-between" amps down
    a00 = qdd_bundle_ptr_amp(QDD_PTR(a00), amp_mul(QDD_AMP(a_low), QDD_AMP(a00)));
    a10 = qdd_bundle_ptr_amp(QDD_PTR(a10), amp_mul(QDD_AMP(a_low), QDD_AMP(a10)));
    a01 = qdd_bundle_ptr_amp(QDD_PTR(a01), amp_mul(QDD_AMP(a_high),QDD_AMP(a01)));
    a11 = qdd_bundle_ptr_amp(QDD_PTR(a11), amp_mul(QDD_AMP(a_high),QDD_AMP(a11)));
    b00 = qdd_bundle_ptr_amp(QDD_PTR(b00), amp_mul(QDD_AMP(b_low), QDD_AMP(b00)));
    b10 = qdd_bundle_ptr_amp(QDD_PTR(b10), amp_mul(QDD_AMP(b_low), QDD_AMP(b10)));
    b01 = qdd_bundle_ptr_amp(QDD_PTR(b01), amp_mul(QDD_AMP(b_high),QDD_AMP(b01)));
    b11 = qdd_bundle_ptr_amp(QDD_PTR(b11), amp_mul(QDD_AMP(b_high),QDD_AMP(b11)));

    // 3. recursive calls (8 tasks: SPAWN 7, CALL 1)
    // |a00 a01| |b00 b01| = b00|a00| + b10|a01| , b01|a00| + b11|a01|
    // |a10 a11| |b10 b11|      |a10|      |a11|      |a10|      |a11|
    QDD a00_b00, a00_b01, a10_b00, a10_b01, a01_b10, a01_b11, a11_b10, a11_b11;
    nextvar++;
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a00, b00, nvars, nextvar)); // 1
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a00, b01, nvars, nextvar)); // 2
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a10, b00, nvars, nextvar)); // 3
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a10, b01, nvars, nextvar)); // 4
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a01, b10, nvars, nextvar)); // 5
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a01, b11, nvars, nextvar)); // 6
    qdd_refs_spawn(SPAWN(qdd_matmat_mult_rec, a11, b10, nvars, nextvar)); // 7
    a11_b11 = CALL(qdd_matmat_mult_rec, a11, b11, nvars, nextvar);
    qdd_refs_push(a11_b11);
    a11_b10 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 7
    a01_b11 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 6
    a01_b10 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 5
    a10_b01 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 4
    a10_b00 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 3
    a00_b01 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 2
    a00_b00 = qdd_refs_sync(SYNC(qdd_matmat_mult_rec)); // 1
    qdd_refs_pop(1);
    nextvar--;

    // 4. gather results of multiplication
    QDD lh1, lh2, rh1, rh2;
    lh1 = qdd_makenode(2*nextvar+1, a00_b00, a10_b00);
    lh2 = qdd_makenode(2*nextvar+1, a01_b10, a11_b10);
    rh1 = qdd_makenode(2*nextvar+1, a00_b01, a10_b01);
    rh2 = qdd_makenode(2*nextvar+1, a01_b11, a11_b11);

    // 5. add resulting qdds
    QDD lh, rh;
    qdd_refs_spawn(SPAWN(qdd_plus, lh1, lh2));
    rh = CALL(qdd_plus, rh1, rh2);
    qdd_refs_push(rh);
    lh = qdd_refs_sync(SYNC(qdd_plus));
    qdd_refs_pop(1);

    // 6. put left and right halves of matix together
    res = qdd_makenode(2*nextvar, lh, rh);

    // Insert in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_MATMAT_MULT, nextvar, QDD_PTR(a), QDD_PTR(b), res)) 
            sylvan_stats_count(QDD_MULT_CACHEDPUT);
    }

    // 7. multiply w/ product of root amps
    AMP prod = amp_mul(QDD_AMP(a), QDD_AMP(b));
    AMP new_root_amp = amp_mul(prod, QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

    return res;
}

QDD
qdd_scalar_mult(QDD qdd, complex_t c)
{
    c = comp_mul(c, comp_value(QDD_AMP(qdd)));
    AMP new_root_amp = qdd_comp_lookup(c);
    return qdd_bundle_ptr_amp(QDD_PTR(qdd), new_root_amp);
}

QDD
qdd_increase_all_vars(QDD qdd, int k)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        return qdd;
    }

    // Check cache
    uint64_t res_ptr;
    if (cache_get3(CACHE_QDD_INC_VARS, QDD_PTR(qdd), k, 0, &res_ptr)) {
        return qdd_bundle_ptr_amp(res_ptr, QDD_AMP(qdd));
    }

    // Get node info
    QDD low, high;
    qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
    qddnode_getchilderen(node, &low, &high);
    BDDVAR curvar = qddnode_getvar(node);

    // Recursive apply to children (TODO: lace?)
    low  = qdd_increase_all_vars(low, k);
    high = qdd_increase_all_vars(high, k);
    QDD res = qdd_makenode(curvar+k, low, high);

    // We assume the input QDDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(QDD_AMP(res) == C_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_QDD_INC_VARS, QDD_PTR(qdd), k, 0, QDD_PTR(res));
    return qdd_bundle_ptr_amp(QDD_PTR(res), QDD_AMP(qdd));
}

QDD
qdd_replace_terminal(QDD a, PTR b)
{
    if (QDD_PTR(a) == QDD_TERMINAL) {
        return qdd_bundle_ptr_amp(b, QDD_AMP(a));
    }

    // Check cache
    uint64_t res_ptr;
    if (cache_get3(CACHE_QDD_REPLACE_TERMINAL, QDD_PTR(a), b, 0, &res_ptr)) {
        return qdd_bundle_ptr_amp(res_ptr, QDD_AMP(a));
    }

    // Get node info
    QDD low, high;
    qddnode_t node = QDD_GETNODE(QDD_PTR(a));
    qddnode_getchilderen(node, &low, &high);
    BDDVAR var = qddnode_getvar(node);

    // Recursive apply to children and makenode (TODO: lace?)
    low  = qdd_replace_terminal(low, b);
    high = qdd_replace_terminal(high, b);
    QDD res = qdd_makenode(var, low, high);
    
    // We assume the input QDDs are already normalized in terms of edge weights
    // so we expect no normalization is needed
    assert(QDD_AMP(res) == C_ONE);

    // Put res (ptr) in cache and return
    cache_put3(CACHE_QDD_REPLACE_TERMINAL, QDD_PTR(a), b, 0, QDD_PTR(res));
    return qdd_bundle_ptr_amp(QDD_PTR(res), QDD_AMP(a));
}

QDD
qdd_tensor_prod(QDD a, QDD b, BDDVAR nvars_a)
{
    // Shift all vars of 'b' by 'nvars_a'
    b = qdd_increase_all_vars(b, nvars_a);
 
    // Stack 'a' on top of 'b' (and multiply root edge of 'a' with that of 'b')
    QDD res = qdd_replace_terminal(a, QDD_PTR(b));
    AMP new_root_amp = amp_mul(QDD_AMP(res), QDD_AMP(b));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

    return res;
}

/******************************</applying gates>*******************************/



/*********************<applying (controlled) sub-circuits>*********************/

QDD
qdd_circuit_swap(QDD qdd, BDDVAR qubit1, BDDVAR qubit2)
{
    assert (qubit1 < qubit2);
    
    LACE_ME; // TODO: make this function a LACE function instead?
    QDD res;

    // CNOT
    res = qdd_cgate(qdd, GATEID_X, qubit1, qubit2);
    // upside down CNOT (equivalent)
    res = qdd_gate(res, GATEID_H, qubit1);
    res = qdd_cgate(res, GATEID_Z, qubit1, qubit2);
    res = qdd_gate(res, GATEID_H, qubit1);
    // CNOT
    res = qdd_cgate(res, GATEID_X, qubit1, qubit2);

    return res;
}

QDD
qdd_circuit_reverse_range(QDD qdd, BDDVAR first, BDDVAR last)
{
    QDD res = qdd;
    BDDVAR a, b;
    int num_qubits = (last - first) + 1;
    for (int j = 0; j < (int)(num_qubits/2); j++) {
        a = first + j;
        b = last  - j;
        res = qdd_circuit_swap(res, a, b);
    }
    return res;
}

QDD
qdd_circuit_QFT(QDD qdd, BDDVAR first, BDDVAR last)
{
    LACE_ME;

    int k;
    QDD res = qdd;
    BDDVAR a, b;
    for (a = first; a <= last; a++) {
        
        // H gate on current qubit
        res = qdd_gate(res, GATEID_H, a);

        // Controlled phase gates on all qubits below
        for (b = a+1; b <= last; b++) {
            k = (b - a) + 1;
            res = qdd_cgate(res, GATEID_Rk(k), a, b);
        }
    }

    // Note that we're not swapping the qubit order in this function

    return res;
}

QDD
qdd_circuit_QFT_inv(QDD qdd, BDDVAR first, BDDVAR last)
{
    LACE_ME;
    
    int k;
    QDD res = qdd;
    BDDVAR a, b;

    // Note that we're not swapping the qubit order in this function
    
    // H gates and phase gates (but now backwards)
    for (a = last + 1; a-- > first; ) { // weird for-loop because BDDVARs are unsigned

        // Controlled phase gates (negative angles this time)
        for (b = last; b >= (a+1); b--){
            k = (b - a) + 1;
            res = qdd_cgate(res, GATEID_Rk_dag(k), a, b);
        }

        // H on current qubit
        res = qdd_gate(res, GATEID_H, a);
    }

    return res;
}

QDD
qdd_circuit(QDD qdd, uint32_t circ_id, BDDVAR t1, BDDVAR t2)
{
    switch (circ_id) {  // don't judge me please
        case CIRCID_swap          : return qdd_circuit_swap(qdd, t1, t2);
        case CIRCID_reverse_range : return qdd_circuit_reverse_range(qdd, t1, t2);
        case CIRCID_QFT           : return qdd_circuit_QFT(qdd, t1, t2);
        case CIRCID_QFT_inv       : return qdd_circuit_QFT_inv(qdd, t1, t2);
        default :
            assert ("Invalid circuit ID" && false);
            return QDD_TERMINAL;
    }
}

TASK_IMPL_6(QDD, qdd_ccircuit, QDD, qdd, uint32_t, circ_id, BDDVAR*, cs, uint32_t, ci, BDDVAR, t1, BDDVAR, t2)
{
    // Cache lookup
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_SUBCIRC, GATE_OPID_40(0, ci, 0), qdd, 
                       GATE_OPID_64(circ_id, cs[0], cs[1], cs[2], t1, t2),
                       &res)) {
            return res;
        }
    }

    // Get current control qubit
    BDDVAR c = cs[ci];

    // If no more control qubits, apply sub-circ here
    if (c == QDD_INVALID_VAR || ci > MAX_CONTROLS) {
        res = qdd_circuit(qdd, circ_id, t1, t2);
        // the gates in qdd_circuit already took care of multiplying the input 
        // root amp with normalization, so no need to do that here again
    }
    else {
        BDDVAR var;
        QDD low, high;
        qdd_get_topvar(qdd, c, &var, &low, &high);
        assert(var <= c);

        if (var == c) {
            ci++; // next control qubit
            high = CALL(qdd_ccircuit, high, circ_id, cs, ci, t1, t2);
            ci--;
        }
        else {
            // recursive call to both children
            qdd_refs_spawn(SPAWN(qdd_ccircuit, high, circ_id, cs, ci, t1, t2));
            low = CALL(qdd_ccircuit, low, circ_id, cs, ci, t1, t2);
            qdd_refs_push(low);
            high = qdd_refs_sync(SYNC(qdd_ccircuit));
            qdd_refs_pop(1);
        }
        res = qdd_makenode(var, low, high);
        // Multiply root amp of sum with input root amp 
        AMP new_root_amp = amp_mul(QDD_AMP(qdd), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    }
    
    // Add to cache, return
    if (cachenow) {
        cache_put3(CACHE_QDD_SUBCIRC, GATE_OPID_40(0, ci, 0), qdd, 
                   GATE_OPID_64(circ_id, cs[0], cs[1], cs[2], t1, t2), 
                   res);
    }
    return res;
}

QDD
qdd_all_control_phase_rec(QDD qdd, BDDVAR k, BDDVAR n, bool *x)
{
    assert(k < n);
    
    bool skipped_k = false;
    qddnode_t node;
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped_k = true;
    }
    else {
        node = QDD_GETNODE(QDD_PTR(qdd));
        BDDVAR var = qddnode_getvar(node);
        if(var > k) {
            skipped_k = true;
        }
    }

    QDD low, high;
    if (skipped_k) {
        // insert skipped node
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
    }
    else {
        // case var == k (var < k shouldn't happen)
        qddnode_getchilderen(node, &low, &high);
    }

    // terminal case, apply phase depending on x[k] (control k on 0 or 1)
    if (k == (n-1)) {
        if (x[k] == 1) {
            AMP new_amp = amp_mul(QDD_AMP(high), C_MIN_ONE);
            high = qdd_bundle_ptr_amp(QDD_PTR(high), new_amp);
        }
        else {
            AMP new_amp = amp_mul(QDD_AMP(low), C_MIN_ONE);
            low = qdd_bundle_ptr_amp(QDD_PTR(low), new_amp);
        }
    }
    // non terminal case, choose low/high depending on x[k] (control k on 0 or 1)
    else {
        if (x[k] == 1) {
            k++; // next level
            high = qdd_all_control_phase_rec(high, k, n, x);
            k--;
        }
        else {
            k++;
            low = qdd_all_control_phase_rec(low, k, n, x);
            k--;
        }
    }

    QDD res = qdd_makenode(k, low, high);

    // multiply by existing edge weight on qdd
    AMP new_root_amp = amp_mul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

QDD
qdd_all_control_phase(QDD qdd, BDDVAR n, bool *x)
{
    return qdd_all_control_phase_rec(qdd, 0, n, x);
}


/********************</applying (controlled) sub-circuits>*********************/



/***********************<measurements and probabilities>***********************/

QDD
qdd_measure_q0(QDD qdd, BDDVAR nvars, int *m, double *p)
{  
    LACE_ME;

    // get probabilities for q0 = |0> and q0 = |1>
    double prob_low, prob_high, prob_root;

    QDD low, high;
    BDDVAR var;
    qdd_get_topvar(qdd, 0, &var, &low, &high);

    if (testing_mode) assert(qdd_is_unitvector(qdd, nvars));

    // TODO: don't use doubles here but allow for mpreal ?
    // (e.g. by using AMPs)
    prob_low  = qdd_unnormed_prob(low,  1, nvars);
    prob_high = qdd_unnormed_prob(high, 1, nvars);
    prob_root = amp_to_prob(QDD_AMP(qdd));
    prob_low  *= prob_root;
    prob_high *= prob_root;
    if (fabs(prob_low + prob_high - 1.0) > 1e-6) {
        printf("WARNING: prob sum = %.10lf (%.5lf + %.5lf)\n", prob_low + prob_high, prob_low, prob_high);
        //assert("probabilities don't sum to 1" && false);
    }

    // flip a coin
    float rnd = ((float)rand())/RAND_MAX;
    *m = (rnd < prob_low) ? 0 : 1;
    *p = prob_low;

    // produce post-measurement state
    AMP norm;
    if (*m == 0) {
        high = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        norm = prob_to_amp(prob_low);
    }
    else {
        low  = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        norm = prob_to_amp(prob_high);
    }

    QDD res = qdd_makenode(0, low, high);

    AMP new_root_amp = amp_mul(QDD_AMP(qdd), QDD_AMP(res));
    new_root_amp     = amp_div(new_root_amp, norm);

    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    res = qdd_remove_global_phase(res);
    return res;
}

QDD
qdd_measure_qubit(QDD qdd, BDDVAR k, BDDVAR nvars, int *m, double *p)
{
    if (k == 0) return qdd_measure_q0(qdd, nvars, m, p);
    qdd = qdd_circuit_swap(qdd, 0, k);
    qdd = qdd_measure_q0(qdd, nvars, m, p);
    qdd = qdd_circuit_swap(qdd, 0, k);
    return qdd;
}

QDD
qdd_measure_all(QDD qdd, BDDVAR n, bool* ms, double *p)
{
    LACE_ME;

    qddnode_t node;
    bool skipped;
    BDDVAR var;
    double prob_low, prob_high, prob_path = 1.0, prob_roots = 1.0;

    for (BDDVAR k=0; k < n; k++) {
        // find relevant node (assuming it should be the next one)
        skipped = false;
        if (QDD_PTR(qdd) == QDD_TERMINAL) {
            skipped = true;
        }
        else {
            node = QDD_GETNODE(QDD_PTR(qdd));
            var = qddnode_getvar(node);
            assert(var >= k);
            if (var > k) skipped = true;
        }
        QDD low, high;
        if (skipped) {
            // if skipped q0 is a don't care, treat separately?
            low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
            high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        }
        else {
            qddnode_getchilderen(node, &low, &high);
        }

        prob_low  = qdd_unnormed_prob(low,  k+1, n);
        prob_high = qdd_unnormed_prob(high, k+1, n);
        prob_roots *= amp_to_prob(QDD_AMP(qdd));
        prob_high = prob_high * prob_roots / prob_path;
        prob_low  = prob_low  * prob_roots / prob_path;

        if (fabs(prob_low + prob_high - 1.0) > amp_store_get_tolerance()) {
            printf("prob sum = %.10lf\n", prob_low + prob_high);
            // printf("probabilities don't sum to 1" && false);
        }

        // flip a coin
        float rnd = ((float)rand())/RAND_MAX;
        ms[k] = (rnd < prob_low) ? 0 : 1;

        // Get next edge
        qdd        = (ms[k] == 0) ? low : high;
        prob_path *= (ms[k] == 0) ? prob_low : prob_high;
    }

    *p = prob_path;
    return qdd_create_basis_state(n, ms);
}

TASK_IMPL_3(double, qdd_unnormed_prob, QDD, qdd, BDDVAR, topvar, BDDVAR, nvars)
{
    assert(topvar <= nvars);

    if (topvar == nvars) {
        assert(QDD_PTR(qdd) == QDD_TERMINAL);
        return amp_to_prob(QDD_AMP(qdd));
    }

    // Look in cache
    bool cachenow = ((topvar % granularity) == 0);
    if (cachenow) {
        uint64_t prob_bits;
        if (cache_get3(CACHE_QDD_PROB, sylvan_false, qdd, QDD_PARAM_PACK_16(topvar, nvars), &prob_bits)) {
            sylvan_stats_count(QDD_PROB_CACHED);
            double_hack_t container = (double_hack_t) prob_bits;
            return container.as_double;
        }
    }

    // Check if the node we want is being skipped
    BDDVAR var;
    QDD low, high;
    qdd_get_topvar(qdd, topvar, &var, &low, &high);

    double prob_low, prob_high, prob_root, prob_res; // "prob" = absolute amps squared
    BDDVAR nextvar = topvar + 1;

    SPAWN(qdd_unnormed_prob, high, nextvar, nvars);
    prob_low  = CALL(qdd_unnormed_prob, low, nextvar, nvars);
    prob_high = SYNC(qdd_unnormed_prob);
    prob_root = amp_to_prob(QDD_AMP(qdd));
    prob_res  = prob_root * (prob_low + prob_high);

    // Put in cache and return
    if (cachenow) {
        double_hack_t container = (double_hack_t) prob_res;
        if (cache_put3(CACHE_QDD_PROB, sylvan_false, qdd, QDD_PARAM_PACK_16(topvar, nvars), container.as_int))
            sylvan_stats_count(QDD_PROB_CACHEDPUT);
    }
    return prob_res;
}

AMP
qdd_get_amplitude(QDD q, bool* basis_state)
{
    AMP res = C_ONE;
    QDD low, high;
    for (;;) {
        res = amp_mul(res, QDD_AMP(q));
        
        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        BDDVAR var     = qddnode_getvar(node);
        qddnode_getchilderen(node, &low, &high);

        // Condition low/high choice on basis state vector[var]
        q = (basis_state[var] == 0) ? low : high;
    }

    return res;
}

complex_t
qdd_get_amplitude_as_complex(QDD q, bool *basis_state)
{
    return comp_value(qdd_get_amplitude(q, basis_state));
}


/***********************<measurements and probabilities>***********************/



/***************************<Initial state creation>***************************/

QDD
qdd_create_all_zero_state(BDDVAR n)
{
    bool x[n];
    for (BDDVAR k=0; k<n; k++) x[k] = 0;
    return qdd_create_basis_state(n, x);
}

QDD
qdd_create_basis_state(BDDVAR n, bool* x)
{
    // start at terminal, and build backwards
    QDD low, high, prev = QDD_TERMINAL;

    for (int k = n-1; k >= 0; k--) {
        if (x[k] == 0) {
            low = qdd_bundle_ptr_amp(QDD_PTR(prev), C_ONE);
            high = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        }
        else if (x[k] == 1) {
            low = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
            high = qdd_bundle_ptr_amp(QDD_PTR(prev), C_ONE);
        }
        // add node to unique table
        prev = qdd_makenode(k, low, high);
    }
    return prev;
}

QDD
qdd_stack_matrix(QDD below, BDDVAR k, gate_id_t gateid)
{
    // This function effectively does a Kronecker product gate \tensor below
    BDDVAR s, t;
    QDD u00, u01, u10, u11, low, high, res;

    // Even + uneven variable are used to encode the 4 values
    s = 2*k;
    t = s + 1;

    // Matrix U = [u00 u01
    //             u10 u11] endoded in a small tree
    u00 = qdd_bundle_ptr_amp(QDD_PTR(below), gates[gateid][0]);
    u10 = qdd_bundle_ptr_amp(QDD_PTR(below), gates[gateid][2]);
    u01 = qdd_bundle_ptr_amp(QDD_PTR(below), gates[gateid][1]);
    u11 = qdd_bundle_ptr_amp(QDD_PTR(below), gates[gateid][3]);
    low  = qdd_makenode(t, u00, u10);
    high = qdd_makenode(t, u01, u11);
    res  = qdd_makenode(s, low, high);

    // Propagate common factor on previous root amp to new root amp
    AMP new_root_amp = amp_mul(QDD_AMP(below), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

QDD
qdd_stack_control(QDD case0, QDD case1, BDDVAR k)
{
    // Effectively does |0><0| \tensor case0 + |1><1| \tensor case1
    BDDVAR s, t;
    QDD u00, u01, u10, u11, low, high, res;

    s = 2*k;
    t = s + 1;

    u00 = case0;
    u10 = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
    u01 = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
    u11 = case1;
    low  = qdd_makenode(t, u00, u10);
    high = qdd_makenode(t, u01, u11);
    res  = qdd_makenode(s, low, high);

    // Weights of case0/case1 already dealt with by qdd_makenode
    return res;
}

QDD
qdd_create_all_identity_matrix(BDDVAR n)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qdd_stack_matrix(prev, k, GATEID_I);
    }
    return prev;
}

QDD
qdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        if ((unsigned int)k == t)
            prev = qdd_stack_matrix(prev, k, gateid);
        else
            prev = qdd_stack_matrix(prev, k, GATEID_I);
    }
    return prev;
}

QDD
qdd_create_single_qubit_gates(BDDVAR n, gate_id_t *gateids)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qdd_stack_matrix(prev, k, gateids[k]);
    }
    return prev;
}

QDD
qdd_create_single_qubit_gates_same(BDDVAR n, gate_id_t gateid)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qdd_stack_matrix(prev, k, gateid);
    }
    return prev;
}

QDD
qdd_create_controlled_gate(BDDVAR n, BDDVAR c, BDDVAR t, gate_id_t gateid)
{
    // for now, assume t > c
    assert(t > c);
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    QDD branch0 = QDD_TERMINAL, branch1 = QDD_TERMINAL;
    for (int k = n-1; k>= 0; k--) {
        if ((unsigned int)k > t || (unsigned int) k < c) {
            prev = qdd_stack_matrix(prev, k, GATEID_I);
        }
        else if ((unsigned int) k == t) {
            branch0 = qdd_stack_matrix(prev, k, GATEID_I);
            branch1 = qdd_stack_matrix(prev, k, gateid);
        }
        else if ((unsigned int) k < t && (unsigned int) k > c) {
            branch0 = qdd_stack_matrix(branch0, k, GATEID_I);
            branch1 = qdd_stack_matrix(branch1, k, GATEID_I);
        }
        else if ((unsigned int) k == c) {
            prev = qdd_stack_control(branch0, branch1, k);
        }
        else {
            assert("all cases should have been covered" && false);
        }
    }
    return prev;
}

QDD
qdd_create_multi_cgate_rec(BDDVAR n, int *c_options, gate_id_t gateid, BDDVAR k)
{
    // (assumes controls above target)
    // c_options[k] = -1 -> ignore qubit k (apply I)
    // c_options[k] =  0 -> control on q_k = |0>
    // c_options[k] =  1 -> control on q_k = |1>
    // c_options[k] =  2 -> target qubit (for now assume 1 target)

    // Terminal case
    if (k == n) {
        return qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    }

    // TODO: catching GATEID_I avoids exp number of recrusive calls, but this 
    // entire function might be handled in a better way(?)
    if (gateid == GATEID_I) {
        QDD identities = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
        for (int j = n-1; j >= (int) k; j--) {
            identities = qdd_stack_matrix(identities, j, GATEID_I);
        }
        return identities;
    }

    // Recursively build matrix
    BDDVAR next_k = k + 1;

    // -1 : Ignore qubit (apply I)
    if (c_options[k] == -1) {
        QDD below = qdd_create_multi_cgate_rec(n, c_options, gateid, next_k);
        return qdd_stack_matrix(below, k, GATEID_I);
    }
    // 0 : control on q_k = |0> (apply gateid to low branch)
    else if (c_options[k] == 0) {
        QDD case0 = qdd_create_multi_cgate_rec(n, c_options, gateid, next_k);
        QDD case1 = qdd_create_multi_cgate_rec(n, c_options, GATEID_I, next_k);
        return qdd_stack_control(case0, case1, k);
    }
    // 1 : control on q_k = |1> (apply gateid to high branch)
    else if (c_options[k] == 1) {
        QDD case0 = qdd_create_multi_cgate_rec(n, c_options, GATEID_I, next_k);
        QDD case1 = qdd_create_multi_cgate_rec(n, c_options, gateid, next_k);
        return qdd_stack_control(case0, case1, k);
    }
    // 2 : target qubit
    else if (c_options[k] == 2) {
        QDD below = qdd_create_multi_cgate_rec(n, c_options, gateid, next_k);
        return qdd_stack_matrix(below, k, gateid);
    }
    else {
        printf("Invalid option %d for qubit %d (options = {-1,0,1,2}\n", c_options[k], k);
        exit(1);
    }
}

QDD
qdd_create_all_control_phase(BDDVAR n, bool *x)
{
    QDD identity = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    QDD ccphase  = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);

    // Start with (-1)Z gate on last qubit. Z if control on 1 and -Z if 0.
    if (x[n-1] == 1) {
        ccphase = qdd_stack_matrix(ccphase, n-1, GATEID_Z);
    }
    else if (x[n-1] == 0) {
        ccphase = qdd_stack_matrix(ccphase, n-1, GATEID_Z);
        ccphase = qdd_bundle_ptr_amp(QDD_PTR(ccphase), amp_neg(QDD_AMP(ccphase)));
    }

    // Stack remaining controls
    for (int k = n-2; k >= 0; k--) {
        // "Identity stack" for doing nothing on each qubit's non-control branch
        identity = qdd_stack_matrix(identity, k+1, GATEID_I);

        // Check if this qubit should be controlled on 0 or 1
        if (x[k] == 1)
            ccphase = qdd_stack_control(identity, ccphase, k);
        else if (x[k] == 0)
            ccphase = qdd_stack_control(ccphase, identity, k);
    }

    return ccphase;
}

/**************************</Initial state creation>***************************/

QDD
qdd_remove_global_phase(QDD qdd)
{
    // Remove global phase by replacing amp of qdd with absolute value of amp
    AMP abs = amp_abs(QDD_AMP(qdd));
    QDD res = qdd_bundle_ptr_amp(QDD_PTR(qdd), abs);
    return res;
}

bool
qdd_equivalent(QDD a, QDD b, int n, bool exact, bool verbose)
{
    bool has_next = true;
    AMP amp_a, amp_b;
    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    while(has_next){
        amp_a = qdd_get_amplitude(a, x);
        amp_b = qdd_get_amplitude(b, x);
        if(exact){
            if(!amp_exact_equal(amp_a, amp_b)){
                if(verbose){
                    // this printing won't work for the mpreal backend
                    // TODO: fix above
                    _print_bitstring(x, n, true);
                    printf(", amp a ="); comp_print(comp_value(amp_a));
                    printf(" != amp b ="); comp_print(comp_value(amp_b));
                    printf("\n");
                }
                return false;
            }
        }
        else{
            if(!amp_approx_equal(amp_a, amp_b)){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", amp a ="); comp_print(comp_value(amp_a));
                    printf(" !~= amp b ="); comp_print(comp_value(amp_b));
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
qdd_is_close_to_unitvector(QDD qdd, BDDVAR n, double tol)
{
    bool WRITE_TO_FILE = false;
    bool has_next = true;
    AMP a;
    bool x[n];
    for(BDDVAR k=0; k<n; k++) x[k] = 0;

    double sum_abs_squares = 0.0;
    while(has_next){
        a = qdd_get_amplitude(qdd, x);
        sum_abs_squares += amp_to_prob(a);
        has_next = _next_bitstring(x, n);
    }

    if (fabs(sum_abs_squares - 1.0) < tol) {
        if (WRITE_TO_FILE) {
            FILE *fp;
            fp = fopen("is_unitvector_true.dot", "w");
            qdd_fprintdot(fp, qdd, false);
            fclose(fp);
        }
        return true;
    }
    else {
        //printf("probs sum to %.30lf\n", sum_abs_squares);
        if (WRITE_TO_FILE) {
            FILE *fp;
            fp = fopen("is_unitvector_false.dot", "w");
            qdd_fprintdot(fp, qdd, false);
            fclose(fp);
        }
        return false;
    }
}

bool
qdd_is_unitvector(QDD qdd, BDDVAR n)
{
    return qdd_is_close_to_unitvector(qdd, n, amp_store_get_tolerance()*10);
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

TASK_IMPL_3(bool, qdd_is_ordered, QDD, qdd, BDDVAR, nextvar, BDDVAR, nvars)
{
    // Terminal case
    if (QDD_PTR(qdd) == QDD_TERMINAL) return true;

    // Get top node
    BDDVAR var;
    QDD low, high;
    qdd_get_topvar(qdd, QDD_INVALID_VAR, &var, &low, &high);

    // Check variable order
    if (var >= nvars || var < nextvar) return false;

    // Check cache
    uint64_t res;
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QDD_IS_ORDERED, QDD_PTR(qdd), sylvan_false, sylvan_false, &res)) {
            return (bool) res;
        }
    }

    // Recursive calls
    bool res_low, res_high;
    var++;
    SPAWN(qdd_is_ordered, high, var, nvars);
    res_low  = CALL(qdd_is_ordered, low, var, nvars);
    res_high = SYNC(qdd_is_ordered);
    var--;
    res = (res_low && res_high);

    // Put res in cache and return
    if (cachenow) cache_put3(CACHE_QDD_IS_ORDERED, QDD_PTR(qdd), sylvan_false, sylvan_false, res);
    return res;
}

/*******************************<logging stats>********************************/

bool qdd_stats_logging = false;
uint32_t statslog_granularity = 1;
uint64_t statslog_buffer = 10; // TODO: remove manual buffer flushing
FILE *qdd_logfile;
uint64_t nodes_peak = 0;
double nodes_avg = 0;
uint64_t logcounter = 0;
uint64_t logtrycounter = 0;

void
qdd_stats_start(FILE *out)
{
    if (out == NULL) return;
    qdd_stats_logging = true;
    qdd_logfile = out;
    fprintf(qdd_logfile, "nodes, amps\n");
    nodes_peak = 0;
    logcounter = 0;
    logtrycounter = 0;
}

void
qdd_stats_set_granularity(uint32_t g)
{
    if (g == 0) statslog_granularity = 1;
    else statslog_granularity = g;
}

void
qdd_stats_log(QDD qdd)
{
    if (!qdd_stats_logging) return;

    // only log every 'statslog_granularity' calls of this function
    if (logtrycounter++ % statslog_granularity != 0) return;

    // Insert info
    uint64_t num_nodes = qdd_countnodes(qdd);
    uint64_t num_amps  = count_amplitude_table_enries();
    fprintf(qdd_logfile, "%ld,%ld\n", num_nodes, num_amps);
    logcounter++;

    // manually flush every 'statslog_buffer' entries
    if (logcounter % statslog_buffer == 0)
        fflush(qdd_logfile);

    // peak nodes
    if (num_nodes > nodes_peak)
        nodes_peak = num_nodes;
    
    // (online) avg nodes
    double a = 1.0/(double)logcounter;
    double b = 1.0 - a;
    nodes_avg = a * (double)num_nodes + b * nodes_avg;
}

uint64_t
qdd_stats_get_nodes_peak()
{
    return nodes_peak;
}

double
qdd_stats_get_nodes_avg()
{
    return nodes_avg;
}

uint64_t
qdd_stats_get_logcounter()
{
    return logtrycounter;
}

void
qdd_stats_finish()
{
    if (!qdd_stats_logging) return;
    fflush(qdd_logfile);
    qdd_stats_logging = false;
    nodes_peak = 0;
    logcounter = 0;
    logtrycounter = 0;
}

/******************************</logging stats>********************************/


/**************************<printing & file writing>***************************/

/**
 * Pretty prints the information contained in `n`.
 */
static void qddnode_pprint(qddnode_t n)
{
    BDDVAR var = qddnode_getvar(n);
    QDD low, high;
    qddnode_getchilderen(n, &low, &high);
    printf("[var=%d, low=%lx, high=%lx, ", 
             var,
             QDD_PTR(low),
             QDD_PTR(high));
    if(QDD_AMP(low) == C_ZERO)      printf("a=C_ZERO, ");
    else if(QDD_AMP(high) == C_ONE)  printf("a=C_ONE, ");
    else {
        printf("a=%lx, ", QDD_AMP(low));
        printf("("); comp_print(comp_value(QDD_AMP(low))); printf(")");
    }                      
    if(QDD_AMP(high) == C_ZERO)     printf("b=C_ZERO ");
    else if(QDD_AMP(high) == C_ONE) printf("b=C_ONE, ");
    else {                     
        printf("b=%lx", QDD_AMP(high));
        printf("("); comp_print(comp_value(QDD_AMP(high))); printf(")");
    }
    printf("]\n");
}

void
_print_qdd(QDD q)
{
    if(QDD_PTR(q) != QDD_TERMINAL){
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        if(!qddnode_getmark(node)){
            qddnode_setmark(node, 1);
            printf("%lx\t", QDD_PTR(q));
            qddnode_pprint(node);
            QDD low, high;
            qddnode_getchilderen(node, &low, &high);
            _print_qdd(low);
            _print_qdd(high);
        }
    }
}

void
qdd_printnodes(QDD q)
{
    printf("root edge: %lx, %lx = ",QDD_PTR(q), QDD_AMP(q));
    comp_print(comp_value(QDD_AMP(q)));
    printf("\n");
    _print_qdd(q);
    qdd_unmark_rec(q);
}

static void
qdd_fprintdot_edge_label(FILE *out, AMP a)
{
    fprintf(out, ", label=\"");
    if (a == C_ONE) {}
    else if (a == C_ZERO) { fprintf(out, "0"); }
    else if (a == C_MIN_ONE) { fprintf(out, "-1"); }
    else {
        complex_t val = comp_value(a);
        if (val.r != 0.0) fprintf(out, "%.2e", (double) val.r);
        if (val.i > 0.0) fprintf(out, "+%.2ei", (double) val.i);
        else if (val.i < 0.0) fprintf(out, "%.2ei", (double) val.i);
    }
    fprintf(out, "\"");
}

static void
qdd_fprintdot_rec(FILE *out, QDD qdd, bool draw_zeros)
{
    // terminal node
    if(QDD_PTR(qdd) == QDD_TERMINAL) return;

    qddnode_t n = QDD_GETNODE(QDD_PTR(qdd));
    if (qddnode_getmark(n)) return;
    qddnode_setmark(n, 1);

    // add this node
    fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n",
            QDD_PTR(qdd), qddnode_getvar(n));

    
    // children of this node
    QDD low, high;
    qddnode_getchilderen(n, &low, &high);
    qdd_fprintdot_rec(out, low, draw_zeros);
    qdd_fprintdot_rec(out, high, draw_zeros);

    // add edge from this node to each child (unless weight 0)
    if (draw_zeros || QDD_AMP(low) != C_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed",
                    QDD_PTR(qdd), QDD_PTR(low));
        qdd_fprintdot_edge_label(out, QDD_AMP(low));
        fprintf(out, "];\n");
    }
    if (draw_zeros || QDD_AMP(high) != C_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                    QDD_PTR(qdd), QDD_PTR(high));
        qdd_fprintdot_edge_label(out, QDD_AMP(high));
        fprintf(out, "];\n");
    }
}

void
qdd_fprintdot(FILE *out, QDD qdd, bool draw_zeros)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir=forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid", QDD_PTR(qdd));
    qdd_fprintdot_edge_label(out, QDD_AMP(qdd));
    fprintf(out, "];\n");

    // terminal node
    fprintf(out, "%lu [shape=box, label=\"T\"];\n", QDD_TERMINAL);

    // recursively add nodes
    qdd_fprintdot_rec(out, qdd, draw_zeros);
    qdd_unmark_rec(qdd);

    fprintf(out, "}\n");
}

/**************************</printing & file writing>**************************/
