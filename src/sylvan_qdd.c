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

//static int granularity = 1; // default
static bool testing_mode = 0; // turns on/off (expensive) sanity checks

/****************< (bit level) manipulation of QDD / qddnode_t >***************/
/**
 * QDD edge structure (64 bits)
 *       1 bit:  marked/unmarked flag (same place as MTBDD)
 *      23 bits: index of edge weight in ctable (AMP)
 *      40 bits: index of next node in node table (PTR)
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
 *      40 bits: low edge pointer to next node (PTR)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag (same place as MTBDD)
 *      23 bits: index of edge weight of high edge in ctable (AMP)
 *      40 bits: high edge pointer to next node (PTR)
 */
static const QDD qdd_marked_mask  = 0x8000000000000000LL;
static const QDD qdd_var_mask_low = 0x7f80000000000000LL;
static const QDD qdd_amp_pos_mask = 0x0040000000000000LL;
static const QDD qdd_amp_val_mask = 0x0020000000000000LL;
static const QDD qdd_amp_mask     = 0x7fffff0000000000LL;
static const QDD qdd_ptr_mask     = 0x000000ffffffffffLL;

typedef struct __attribute__((packed)) qddnode {
    QDD low, high;
} *qddnode_t; // 16 bytes


/**
 * Gets only the AMP information of a QDD edge `q`.
 */
static inline AMP
QDD_AMP(QDD q)
{
    return (q & qdd_amp_mask) >> 40; // 23 bits
}

/**
 * Gets only the PTR information of a QDD edge `q`.
 */
static inline PTR
QDD_PTR(QDD q)
{
    return q & qdd_ptr_mask; // 40 bits
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
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= (1<<23));
    return (a << 40 | p);
}

static void
qddnode_getchild_ptrs_amps(qddnode_t n, PTR *low, PTR *high, AMP *a, AMP *b)
{
    *low  = qddnode_getptrlow(n);
    *high = qddnode_getptrhigh(n);

    bool norm_pos = (n->low & qdd_amp_pos_mask) >> 54;
    bool norm_val = (n->low & qdd_amp_val_mask) >> 53;
    if (norm_pos == 0) { // low amp is C_ZERO or C_ONE, high amp in ctable
        *a = (norm_val == 0) ? C_ZERO : C_ONE;
        *b = QDD_AMP(n->high);
    }
    else { // high amp is C_ZERO or C_ONE, low amp in ctable
        *b = (norm_val == 0) ? C_ZERO : C_ONE;
        *a = QDD_AMP(n->high);
    }
}

static void
qddnode_getchilderen(qddnode_t n, QDD *low, QDD *high)
{
    PTR l, h;
    AMP a, b;
    qddnode_getchild_ptrs_amps(n, &l, &h, &a, &b);
    *low  = qdd_bundle_ptr_amp(l, a);
    *high = qdd_bundle_ptr_amp(h, b);
}

static void
qddnode_pack(qddnode_t n, BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    assert(a == C_ZERO || a == C_ONE || b == C_ZERO || b == C_ONE);

    AMP amp_high;
    bool norm_pos = (a == C_ZERO || a == C_ONE) ? 0 : 1;
    bool norm_val;
    if (norm_pos == 0) {
        norm_val = (a == C_ZERO) ? 0 : 1;
        amp_high = b;
    }
    else {
        norm_val = (b == C_ZERO) ? 0 : 1;
        amp_high = a;
    }

    n->low  = ((uint64_t)var)<<55 | ((uint64_t)norm_pos)<<54 | ((uint64_t)norm_val)<<53 | low;
    n->high = amp_high<<40 | high;
}

// Container for disguising doubles as ints so they can go in Sylvan's cache
// (see also union "hack" in mtbdd_satcount)
typedef union {
    double   as_double;
    uint64_t as_int;
} double_hack_t;

typedef union {
    complex_t as_comp;
    uint64_t  as_int[2];
} comp_hack_t;

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
sylvan_init_qdd(size_t ctable_size)
{
    if (qdd_initialized) return;
    qdd_initialized = 1;

    sylvan_register_quit(qdd_quit);
    sylvan_gc_add_mark(TASK(qdd_gc_mark_external_refs));
    sylvan_gc_add_mark(TASK(qdd_gc_mark_protected));

    refs_create(&qdd_refs, 1024);
    if (!qdd_protected_created) {
        protect_create(&qdd_protected, 4096);
        qdd_protected_created = 1;
    }

    init_amplitude_table(ctable_size);

    LACE_ME;
    CALL(qdd_refs_init);
}

void
qdd_set_testing_mode(bool on)
{
    testing_mode = on;
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

static AMP __attribute__((unused))
qdd_amp_normalize_low(AMP *low, AMP *high)
{
    // Normalize using low if low != 0
    AMP norm;
    if(*low != C_ZERO){
        complex_t cl = comp_value(*low);
        complex_t ch = comp_value(*high);
        ch    = comp_div(ch, cl);
        *high = qdd_comp_lookup(ch);
        norm  = *low;
        *low  = C_ONE;
    }
    else {
        norm  = *high;
        *high = C_ONE;
    }
    return norm;
}

static AMP __attribute__((unused))
qdd_amp_normalize_largest(AMP *low, AMP *high)
{
    AMP norm;
    if (*low == *high) {
        norm  = *low;
        *low  = C_ONE;
        *high = C_ONE;
        return norm;
    }

    // Normalize using the absolute greatest value
    complex_t cl = comp_value(*low);
    complex_t ch = comp_value(*high);
    if ( (cl.r*cl.r + cl.i*cl.i)  >=  (ch.r*ch.r + ch.i*ch.i) ) {
        ch = comp_div(ch, cl);
        *high = qdd_comp_lookup(ch);
        norm = *low;
        *low  = C_ONE;
    }
    else {
        cl = comp_div(cl, ch);
        *low = qdd_comp_lookup(cl);
        norm  = *high;
        *high = C_ONE;
    }
    return norm;
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
        AMP norm = qdd_amp_normalize_largest(&low_amp, &high_amp);
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

static int auto_gc_ctable     = 1;
static double ctable_gc_thres = 0.5;

void
qdd_set_auto_gc_ctable(bool enabled)
{
    auto_gc_ctable = enabled;
}

void
qdd_set_gc_ctable_thres(double fraction_filled)
{
    ctable_gc_thres = fraction_filled;
}


void
qdd_gc_ctable(QDD *keep)
{
    LACE_ME;
    // 1. Create new amp table
    init_new_empty_table();

    // 2. Fill new table with amps present in given QDDs
    //for (int i = 0; i < n_qdds; i++) qdds[i] = _fill_new_amp_table(qdds[i]);
    *keep = _fill_new_amp_table(*keep);

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

void
qdd_test_gc_ctable(QDD *keep)
{
    uint64_t entries = get_ctable_entries_estimate();
    uint64_t size    = get_ctable_size();
    if ( ((double)entries / (double)size) > ctable_gc_thres )
        qdd_gc_ctable(keep);
}

/*************************</cleaning amplitude table>**************************/



/*******************************<applying gates>*******************************/

static void
qdd_do_before_gate(QDD* qdd)
{
    // check if ctable needs gc
    if (auto_gc_ctable) qdd_test_gc_ctable(qdd);

    // log stuff (if logging is enabled)
    qdd_stats_log(*qdd);
}

/* Wrapper for applying a single qubit gate. */
TASK_IMPL_3(QDD, qdd_gate, QDD, qdd, uint32_t, gate, BDDVAR, target)
{
    qdd_do_before_gate(&qdd);
    return qdd_gate_rec(qdd, gate, target);
}

/* Wrapper for applying a controled gate with 1 control qubit. */
TASK_IMPL_4(QDD, qdd_cgate, QDD, qdd, uint32_t, gate, BDDVAR, c, BDDVAR, t)
{
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c, QDD_INVALID_VAR, QDD_INVALID_VAR, QDD_INVALID_VAR};
    return qdd_cgate_rec(qdd, gate, cs, t);
}

/* Wrapper for applying a controled gate with 2 control qubits. */
TASK_IMPL_5(QDD, qdd_cgate2, QDD, qdd, uint32_t, gate, BDDVAR, c1, BDDVAR, c2, BDDVAR, t)
{
    assert(false);
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c1, c2, QDD_INVALID_VAR, QDD_INVALID_VAR};
    return qdd_cgate_rec(qdd, gate, cs, t);
}

/* Wrapper for applying a controled gate with 3 control qubits. */
TASK_IMPL_6(QDD, qdd_cgate3, QDD, qdd, uint32_t, gate, BDDVAR, c1, BDDVAR, c2, BDDVAR, c3, BDDVAR, t)
{
    assert(false);
    qdd_do_before_gate(&qdd);
    BDDVAR cs[4] = {c1, c2, c3, QDD_INVALID_VAR}; // last pos is a buffer
    return qdd_cgate_rec(qdd, gate, cs, t);
}

TASK_IMPL_2(QDD, qdd_plus_amp, QDD, a, QDD, b)
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
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_PLUS, sylvan_false, a, b, &res)) {
            sylvan_stats_count(QDD_PLUS_CACHED);
            return res;
        }
    }

    // If not base/terminal case, pass edge weight of current edge down
    AMP amp_la, amp_ha, amp_lb, amp_hb;
    amp_la = amp_mul(QDD_AMP(a), QDD_AMP(low_a));
    amp_ha = amp_mul(QDD_AMP(a), QDD_AMP(high_a));
    amp_lb = amp_mul(QDD_AMP(b), QDD_AMP(low_b));
    amp_hb = amp_mul(QDD_AMP(b), QDD_AMP(high_b));
    low_a  = qdd_bundle_ptr_amp(QDD_PTR(low_a),  amp_la);
    high_a = qdd_bundle_ptr_amp(QDD_PTR(high_a), amp_ha);
    low_b  = qdd_bundle_ptr_amp(QDD_PTR(low_b),  amp_lb);
    high_b = qdd_bundle_ptr_amp(QDD_PTR(high_b), amp_hb);

    // Recursive calls down
    QDD low, high;
    qdd_refs_spawn(SPAWN(qdd_plus_amp, high_a, high_b));
    low = CALL(qdd_plus_amp, low_a, low_b);
    qdd_refs_push(low);
    high = qdd_refs_sync(SYNC(qdd_plus_amp));
    qdd_refs_pop(1);

    // Put in cache, return
    res = qdd_makenode(topvar, low, high);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_PLUS, sylvan_false, a, b, res)) 
            sylvan_stats_count(QDD_PLUS_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_2(QDD, qdd_plus_comp_wrap, QDD, a, QDD, b)
{
    complex_t ca, cb;
    ca = comp_value(QDD_AMP(a));
    cb = comp_value(QDD_AMP(b));
    return CALL(qdd_plus_complex, QDD_PTR(a), QDD_PTR(b), ca, cb);
}

/**
 * Version of qdd_plus which propagates complex_t's down instead of first
 * hashing them into AMPs and now only hash them into the complex table when
 * nodes are created. The goal of this is to alleviate access to the complex 
 * table (which is shared between threads). Drawback is that now these complex 
 * structs (made of floating point values) are keys for the operation cache.
 */
TASK_IMPL_4(QDD, qdd_plus_complex, PTR, a, PTR, b, complex_t, ca, complex_t, cb)
{
    // Trivial cases // Q: use exact or approx?
    if (comp_exact_equal(ca, comp_zero())) 
        return qdd_bundle_ptr_amp(b, qdd_comp_lookup(cb));
    if (comp_exact_equal(cb, comp_zero()))
        return qdd_bundle_ptr_amp(a, qdd_comp_lookup(ca));

    // Get var(a) and var(b)
    QDD low_a, low_b, high_a, high_b, res;
    BDDVAR var_a = UINT32_MAX, var_b = UINT32_MAX, topvar;
    if (a != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(a);
        var_a  = qddnode_getvar(node);
    }
    if (b != QDD_TERMINAL) {
        qddnode_t node = QDD_GETNODE(b);
        var_b  = qddnode_getvar(node);
    }

    // For both a and b, get children of node with var=top{topvar(a),topvar(b)}
    // TODO: maybe make topvar function which doesn't bundle PTR and AMP, atm 
    // unnecessary bundling/unbundling is happening
    qdd_get_topvar(a, var_b, &topvar, &low_a, &high_a);
    qdd_get_topvar(b, var_a, &topvar, &low_b, &high_b);

    // Base/terminal case: same target and same variable
    if(a == b && var_a == var_b){
        AMP sum = qdd_comp_lookup(comp_add(ca, cb));
        res = qdd_bundle_ptr_amp(a, sum);
        return res;
    }

    // Check cache
    bool cachenow = 1;
    if (cachenow) {
        QDD blank;
        comp_hack_t hca = (comp_hack_t) ca;
        comp_hack_t hcb = (comp_hack_t) cb;
        if (cache_get6((CACHE_QDD_PLUS | a), hca.as_int[0], hca.as_int[1],
                        b, hcb.as_int[0], hcb.as_int[1],
                        &res, &blank)) {
            sylvan_stats_count(QDD_PLUS_CACHED);
            return res;
        }
    }

    // If not base/terminal case, pass edge weight of current edge down
    complex_t comp_la, comp_ha, comp_lb, comp_hb;
    comp_la = comp_mul(ca, comp_value(QDD_AMP(low_a)));
    comp_ha = comp_mul(ca, comp_value(QDD_AMP(high_a)));
    comp_lb = comp_mul(cb, comp_value(QDD_AMP(low_b)));
    comp_hb = comp_mul(cb, comp_value(QDD_AMP(high_b)));

    // Recursive calls down
    QDD low, high;
    qdd_refs_spawn(SPAWN(qdd_plus_complex, QDD_PTR(high_a), QDD_PTR(high_b), comp_ha, comp_hb));
    low = CALL(qdd_plus_complex, QDD_PTR(low_a), QDD_PTR(low_b), comp_la, comp_lb);
    qdd_refs_push(low);
    high = qdd_refs_sync(SYNC(qdd_plus_complex));
    qdd_refs_pop(1);

    // Put in cache, return
    res = qdd_makenode(topvar, low, high);
    if (cachenow) {
        comp_hack_t hca = (comp_hack_t) ca;
        comp_hack_t hcb = (comp_hack_t) cb;
        cache_put6((CACHE_QDD_PLUS | a), hca.as_int[0], hca.as_int[1],
                    b, hcb.as_int[0], hcb.as_int[1],
                    res, 0LL);
    }
    return res;
}

TASK_IMPL_3(QDD, qdd_gate_rec_amp, QDD, q, uint32_t, gate, BDDVAR, target)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_GATE, GATE_OPID_40(gate, target, 0), q, sylvan_false, &res)) {
            sylvan_stats_count(QDD_GATE_CACHED);
            return res;
        }
    }

    BDDVAR var;
    QDD low, high;
    qdd_get_topvar(q, target, &var, &low, &high);
    assert(var <= target);

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
        qdd_refs_spawn(SPAWN(qdd_plus_amp, high1, high2));
        low = CALL(qdd_plus_amp, low1, low2);
        qdd_refs_push(low);
        high = SYNC(qdd_plus_amp);
        res = qdd_makenode(target, low, high);
       
    }
    else { // var < target: not at target qubit yet, recursive calls down
        qdd_refs_spawn(SPAWN(qdd_gate_rec_amp, high, gate, target));
        low = CALL(qdd_gate_rec_amp, low, gate, target);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_gate_rec_amp));
        qdd_refs_pop(1);
        res  = qdd_makenode(var, low, high);
    }

    // Multiply root amp of sum with input root amp, add to cache, return
    AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_GATE, GATE_OPID_40(gate, target, 0), q, sylvan_false, res)) 
            sylvan_stats_count(QDD_GATE_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_3(QDD, qdd_gate_rec_complex, QDD, qdd, uint32_t, gateid, BDDVAR, target)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_GATE, GATE_OPID_40(gateid, target, 0), qdd, sylvan_false, &res)) {
            sylvan_stats_count(QDD_GATE_CACHED);
            return res;
        }
    }

    BDDVAR var;
    QDD low, high;
    qdd_get_topvar(qdd, target, &var, &low, &high);
    assert(var <= target);
    complex_t ca = comp_value(QDD_AMP(low));
    complex_t cb = comp_value(QDD_AMP(high));

    if (var == target) {
        complex_t a_u00, a_u10, b_u01, b_u11;
        a_u00 = comp_mul(ca, comp_value(gates[gateid][0])); // TODO: keep array
        a_u10 = comp_mul(ca, comp_value(gates[gateid][2])); // of gate values
        b_u01 = comp_mul(cb, comp_value(gates[gateid][1])); // as complex_t's ?
        b_u11 = comp_mul(cb, comp_value(gates[gateid][3]));
        QDD res_low, res_high;
        qdd_refs_spawn(SPAWN(qdd_plus_complex, QDD_PTR(low), QDD_PTR(high), a_u10, b_u11));
        res_low = CALL(qdd_plus_complex, QDD_PTR(low), QDD_PTR(high), a_u00, b_u01);
        qdd_refs_push(low);
        res_high = SYNC(qdd_plus_complex);
        res = qdd_makenode(target, res_low, res_high);
    }
    else { // var < target: not at target qubit yet, recursive calls down
        qdd_refs_spawn(SPAWN(qdd_gate_rec_complex, high, gateid, target));
        low = CALL(qdd_gate_rec_complex, low, gateid, target);
        qdd_refs_push(low);
        high = qdd_refs_sync(SYNC(qdd_gate_rec_complex));
        qdd_refs_pop(1);
        res  = qdd_makenode(var, low, high);
    }

    // Multiply root amp of sum with input root amp, add to cache, return
    complex_t c = comp_mul(comp_value(QDD_AMP(qdd)), comp_value(QDD_AMP(res)));
    AMP new_root_amp = qdd_comp_lookup(c);
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    if (cachenow) {
        if (cache_put3(CACHE_QDD_GATE, GATE_OPID_40(gateid, target, 0), qdd, sylvan_false, res))
            sylvan_stats_count(QDD_GATE_CACHEDPUT);
    }
    return res;
}

TASK_IMPL_5(QDD, qdd_cgate_rec_amp, QDD, q, uint32_t, gate, BDDVAR*, cs, uint32_t, ci, BDDVAR, t)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CGATE, sylvan_false, q,
                       GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                       &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            return res;
        }
    }

    // Get current control qubit
    BDDVAR c = cs[ci];

    // If no more control qubits, apply gate here
    if (c == QDD_INVALID_VAR || ci > MAX_CONTROLS) {
        res = CALL(qdd_gate_rec_amp, q, gate, t);
    }
    else {
        BDDVAR var;
        QDD low, high;
        qdd_get_topvar(q, c, &var, &low, &high);
        assert(var <= c);

        // If current node is (one of) the control qubit(s), 
        // control on q_c = |1> (high edge)
        if (var == c) {
            ci++;
            high = CALL(qdd_cgate_rec_amp, high, gate, cs, ci, t);
            ci--;
        }
        // Not at control qubit yet, apply to both childeren.
        else {
            qdd_refs_spawn(SPAWN(qdd_cgate_rec_amp, high, gate, cs, ci, t));
            low = CALL(qdd_cgate_rec_amp, low, gate, cs, ci, t);
            qdd_refs_push(low);
            high = qdd_refs_sync(SYNC(qdd_cgate_rec_amp));
            qdd_refs_pop(1);
        }
        res = qdd_makenode(var, low, high);

        // Multiply root amp of sum with input root amp, add to cache, return
        AMP new_root_amp = amp_mul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    }

    // Add to cache, return
    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE, sylvan_false, q,
                       GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                       res)) {
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
        }
    }
    return res;
}

TASK_IMPL_5(QDD, qdd_cgate_rec_complex, QDD, q, uint32_t, gate, BDDVAR*, cs, uint32_t, ci, BDDVAR, t)
{
    // Check cache
    QDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_CGATE, sylvan_false, q,
                       GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                       &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            return res;
        }
    }

    // Get current control qubit
    BDDVAR c = cs[ci];

    // If no more control qubits, apply gate here
    if (c == QDD_INVALID_VAR || ci > MAX_CONTROLS) {
        res = CALL(qdd_gate_rec_complex, q, gate, t);
    }
    else {
        BDDVAR var;
        QDD low, high;
        qdd_get_topvar(q, c, &var, &low, &high);
        assert(var <= c);

        // If current node is (one of) the control qubit(s), 
        // control on q_c = |1> (high edge)
        if (var == c) {
            ci++;
            high = CALL(qdd_gate_rec_complex, high, gate, t);
            ci--;
        }
        // Not at control qubit yet, apply to both childeren.
        else {
            qdd_refs_spawn(SPAWN(qdd_cgate_rec_complex, high, gate, cs, ci, t));
            low = CALL(qdd_cgate_rec_complex, low, gate, cs, ci, t);
            qdd_refs_push(low);
            high = qdd_refs_sync(SYNC(qdd_cgate_rec_complex));
            qdd_refs_pop(1);
        }
        res = qdd_makenode(var, low, high);

        // Multiply root amp of sum with input root amp, add to cache, return
        complex_t comp = comp_mul(comp_value(QDD_AMP(q)), comp_value(QDD_AMP(res)));
        AMP new_root_amp = qdd_comp_lookup(comp);
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    }

    // Add to cache, return
    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE, sylvan_false, q,
                       GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                       res)) {
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
        }
    }
    return res;
}

TASK_IMPL_4(QDD, qdd_matvec_mult, QDD, mat, QDD, vec, BDDVAR, nvars, BDDVAR, nextvar)
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
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_MATVEC_MULT, sylvan_false, mat, vec, &res)) {
            sylvan_stats_count(QDD_MULT_CACHED);
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

    // 2. pass weights of current edges down
    AMP amp_vec_low, amp_vec_high, amp_u00, amp_u10, amp_u01, amp_u11;
    amp_vec_low  = amp_mul(QDD_AMP(vec), QDD_AMP(vec_low));
    amp_vec_high = amp_mul(QDD_AMP(vec), QDD_AMP(vec_high));
    amp_u00      = amp_mul(amp_mul(QDD_AMP(mat), QDD_AMP(mat_low)), QDD_AMP(u00));
    amp_u10      = amp_mul(amp_mul(QDD_AMP(mat), QDD_AMP(mat_low)), QDD_AMP(u10));
    amp_u01      = amp_mul(amp_mul(QDD_AMP(mat), QDD_AMP(mat_high)),QDD_AMP(u01));
    amp_u11      = amp_mul(amp_mul(QDD_AMP(mat), QDD_AMP(mat_high)),QDD_AMP(u11));
    vec_low  = qdd_bundle_ptr_amp(QDD_PTR(vec_low), amp_vec_low);
    vec_high = qdd_bundle_ptr_amp(QDD_PTR(vec_high),amp_vec_high);
    u00      = qdd_bundle_ptr_amp(QDD_PTR(u00), amp_u00);
    u10      = qdd_bundle_ptr_amp(QDD_PTR(u10), amp_u10);
    u01      = qdd_bundle_ptr_amp(QDD_PTR(u01), amp_u01);
    u11      = qdd_bundle_ptr_amp(QDD_PTR(u11), amp_u11);

    // 3. recursive calls (4 tasks: SPAWN 3, CALL 1)
    // |u00 u01| |vec_low | = vec_low|u00| + vec_high|u01|
    // |u10 u11| |vec_high|          |u10|           |u11|
    QDD res_low00, res_low10, res_high01, res_high11;
    nextvar++;
    qdd_refs_spawn(SPAWN(qdd_matvec_mult, u00, vec_low,  nvars, nextvar)); // 1
    qdd_refs_spawn(SPAWN(qdd_matvec_mult, u10, vec_low,  nvars, nextvar)); // 2
    qdd_refs_spawn(SPAWN(qdd_matvec_mult, u01, vec_high, nvars, nextvar)); // 3
    res_high11 = CALL(qdd_matvec_mult, u11, vec_high, nvars, nextvar);
    qdd_refs_push(res_high11);
    res_high01 = qdd_refs_sync(SYNC(qdd_matvec_mult)); // 3
    res_low10  = qdd_refs_sync(SYNC(qdd_matvec_mult)); // 2
    res_low00  = qdd_refs_sync(SYNC(qdd_matvec_mult)); // 1
    qdd_refs_pop(1);
    nextvar--;

    // 4. gather results of multiplication
    QDD res_low, res_high;
    res_low  = qdd_makenode(nextvar, res_low00,  res_low10);
    res_high = qdd_makenode(nextvar, res_high01, res_high11);

    // 5. add resulting qdds
    res = CALL(qdd_plus_amp, res_low, res_high);

    // Insert in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_MATVEC_MULT, sylvan_false, mat, vec, res)) 
            sylvan_stats_count(QDD_MULT_CACHEDPUT);
    }

    return res;
}

TASK_IMPL_4(QDD, qdd_matmat_mult, QDD, a, QDD, b, BDDVAR, nvars, BDDVAR, nextvar)
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
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QDD_MATMAT_MULT, sylvan_false, a, b, &res)) {
            sylvan_stats_count(QDD_MULT_CACHED);
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

    // 2. pass weights of current edges down TODO: use loop
    AMP amp_a00, amp_a10, amp_a01, amp_a11, amp_b00, amp_b10, amp_b01, amp_b11;
    amp_a00 = amp_mul(amp_mul(QDD_AMP(a), QDD_AMP(a_low)), QDD_AMP(a00));
    amp_a10 = amp_mul(amp_mul(QDD_AMP(a), QDD_AMP(a_low)), QDD_AMP(a10));
    amp_a01 = amp_mul(amp_mul(QDD_AMP(a), QDD_AMP(a_high)),QDD_AMP(a01));
    amp_a11 = amp_mul(amp_mul(QDD_AMP(a), QDD_AMP(a_high)),QDD_AMP(a11));
    amp_b00 = amp_mul(amp_mul(QDD_AMP(b), QDD_AMP(b_low)), QDD_AMP(b00));
    amp_b10 = amp_mul(amp_mul(QDD_AMP(b), QDD_AMP(b_low)), QDD_AMP(b10));
    amp_b01 = amp_mul(amp_mul(QDD_AMP(b), QDD_AMP(b_high)),QDD_AMP(b01));
    amp_b11 = amp_mul(amp_mul(QDD_AMP(b), QDD_AMP(b_high)),QDD_AMP(b11));
    a00 = qdd_bundle_ptr_amp(QDD_PTR(a00), amp_a00);
    a10 = qdd_bundle_ptr_amp(QDD_PTR(a10), amp_a10);
    a01 = qdd_bundle_ptr_amp(QDD_PTR(a01), amp_a01);
    a11 = qdd_bundle_ptr_amp(QDD_PTR(a11), amp_a11);
    b00 = qdd_bundle_ptr_amp(QDD_PTR(b00), amp_b00);
    b10 = qdd_bundle_ptr_amp(QDD_PTR(b10), amp_b10);
    b01 = qdd_bundle_ptr_amp(QDD_PTR(b01), amp_b01);
    b11 = qdd_bundle_ptr_amp(QDD_PTR(b11), amp_b11);

    // 3. recursive calls (8 tasks: SPAWN 7, CALL 1)
    // |a00 a01| |b00 b01| = b00|a00| + b10|a01| , b01|a00| + b11|a01|
    // |a10 a11| |b10 b11|      |a10|      |a11|      |a10|      |a11|
    QDD a00_b00, a00_b01, a10_b00, a10_b01, a01_b10, a01_b11, a11_b10, a11_b11;
    nextvar++;
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a00, b00, nvars, nextvar)); // 1
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a00, b01, nvars, nextvar)); // 2
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a10, b00, nvars, nextvar)); // 3
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a10, b01, nvars, nextvar)); // 4
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a01, b10, nvars, nextvar)); // 5
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a01, b11, nvars, nextvar)); // 6
    qdd_refs_spawn(SPAWN(qdd_matmat_mult, a11, b10, nvars, nextvar)); // 7
    a11_b11 = CALL(qdd_matmat_mult, a11, b11, nvars, nextvar);
    qdd_refs_push(a11_b11);
    a11_b10 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 7
    a01_b11 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 6
    a01_b10 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 5
    a10_b01 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 4
    a10_b00 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 3
    a00_b01 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 2
    a00_b00 = qdd_refs_sync(SYNC(qdd_matmat_mult)); // 1
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
    qdd_refs_spawn(SPAWN(qdd_plus_amp, lh1, lh2));
    rh = CALL(qdd_plus_amp, rh1, rh2);
    qdd_refs_push(rh);
    lh = qdd_refs_sync(SYNC(qdd_plus_amp));
    qdd_refs_pop(1);

    // 6. put left and right halves of matix together
    res = qdd_makenode(2*nextvar, lh, rh);

    // Insert in cache
    if (cachenow) {
        if (cache_put3(CACHE_QDD_MATMAT_MULT, sylvan_false, a, b, res)) 
            sylvan_stats_count(QDD_MULT_CACHEDPUT);
    }

    return res;
}

QDD
qdd_scalar_mult(QDD qdd, complex_t c)
{
    c = comp_mul(c, comp_value(QDD_AMP(qdd)));
    AMP new_root_amp = qdd_comp_lookup(c);
    return qdd_bundle_ptr_amp(QDD_PTR(qdd), new_root_amp);
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
    for (a = last + 1; a-- > first; ) { // weird for loop because BDDVARs are unsigned

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
    // TODO? caching here?
    switch (circ_id) {  // don't judge me please
        case CIRCID_swap          : return qdd_circuit_swap(qdd, t1, t2);
        case CIRCID_reverse_range : return qdd_circuit_reverse_range(qdd, t1, t2);
        case CIRCID_QFT           : return qdd_circuit_QFT(qdd, t1, t2);
        case CIRCID_QFT_inv       : return qdd_circuit_QFT_inv(qdd, t1, t2);
        case CIRCID_phi_add_a     : return qdd_phi_add(qdd, t1, t2, shor_bits_a);
        case CIRCID_phi_add_N     : return qdd_phi_add(qdd, t1, t2, shor_bits_N);
        case CIRCID_phi_add_a_inv : return qdd_phi_add_inv(qdd, t1, t2, shor_bits_a);
        case CIRCID_phi_add_N_inv : return qdd_phi_add_inv(qdd, t1, t2, shor_bits_N);
        default :
            assert ("Invalid circuit ID" && false);
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
        // TODO: (I guess this needs to be done every time after makenode, 
        // so maybe put this functionality in makenode function?)
        complex_t comp = comp_mul(comp_value(QDD_AMP(qdd)), comp_value(QDD_AMP(res)));
        AMP new_root_amp = qdd_comp_lookup(comp);
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
            complex_t c = comp_mul(comp_value(QDD_AMP(high)), comp_minus_one());
            AMP new_amp = qdd_comp_lookup(c);
            high = qdd_bundle_ptr_amp(QDD_PTR(high), new_amp);
        }
        else {
            complex_t c = comp_mul(comp_value(QDD_AMP(low)), comp_minus_one());
            AMP new_amp = qdd_comp_lookup(c);
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
    complex_t c = comp_mul(comp_value(QDD_AMP(qdd)), comp_value(QDD_AMP(res)));
    AMP new_root_amp = qdd_comp_lookup(c);
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}

QDD
qdd_all_control_phase(QDD qdd, BDDVAR n, bool *x)
{
    return qdd_all_control_phase_rec(qdd, 0, n, x);
}


/********************</applying (controlled) sub-circuits>*********************/



/*******************************<Shor components>******************************/

uint64_t 
inverse_mod(uint64_t a, uint64_t N) {
    int t = 0;
    int newt = 1;
    int r = N;
    int newr = a;
    int h;
    while(newr != 0) {
        int quotient = r / newr;
        h = t;
        t = newt;
        newt = h - quotient * newt;
        h = r;
        r = newr;
        newr = h - quotient * newr;
    }
    if(r > 1)
        printf("ERROR: a is not invertible\n");
    if(t < 0)
        t = t + N;
    return t;
}

QDD 
qdd_phi_add(QDD qdd, BDDVAR first, BDDVAR last, bool* a) 
{
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = 0; j < (num_qubits-i); j++) {
            if (a[j] == 1) {
                k = num_qubits-j-i; // Rk = 2pi/2^k rotation
                res = qdd_gate(res, GATEID_Rk(k), qubit);
            }
        }
    }
    return res;
}

QDD 
qdd_phi_add_inv(QDD qdd, BDDVAR first, BDDVAR last, bool* a) 
{
    // These are all phase gates, so they'll all commute, so this is the exact
    // same function als qdd_phi_add() but with inversed angles.
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = 0; j < (num_qubits-i); j++) {
            if (a[j] == 1) {
                k = num_qubits-j-i; // Rk_dag = -2pi/2^k rotation
                res = qdd_gate(res, GATEID_Rk_dag(k), qubit);
            }
        }
    }
    return res;
}

QDD
qdd_phi_add_mod(QDD qdd, BDDVAR* cs, uint64_t a, uint64_t N)
{
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)

    // 1.  controlled(c1,c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 2.  phi_add_inv(N)
    qdd = qdd_circuit(qdd, CIRCID_phi_add_N_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 3.  QFT_inv
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 4.  CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 5.  QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 6.  controlled phi_add(N) (control = helper)
    BDDVAR controls[] = {shor_wires.helper, QDD_INVALID_VAR, QDD_INVALID_VAR};
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_N, controls, shor_wires.targ_first, shor_wires.targ_last);
    // 7. controlled(c1, c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 8.  QFT_inv
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 9.  X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 11. X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 12. QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 13. controlled(c1,c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_phi_add_mod_inv(QDD qdd, BDDVAR* cs, uint64_t a, uint64_t N)
{
    // Inverse of function above
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)

    // 13. controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 12. QFT^-1
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 11. X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT^-1 = CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 9.  X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 8.  (QFT^-1)^-1 = QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 7.  controlled(c1, c2) phi_add(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a, cs, shor_wires.targ_first, shor_wires.targ_last);
    // 6.  controlled phi_add_inv(N) (control = helper)
    BDDVAR controls[] = {shor_wires.helper, QDD_INVALID_VAR, QDD_INVALID_VAR};
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_N_inv, controls, shor_wires.targ_first, shor_wires.targ_last);
    // 5.  QFT^-1
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);
    // 4. CNOT^-1 = CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 3. (QFT^-1)^-1 = QFT
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    // 2.  phi_add(N)
    qdd = qdd_circuit(qdd, CIRCID_phi_add_N, shor_wires.targ_first, shor_wires.targ_last);
    // 1.  controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_ccircuit(qdd, CIRCID_phi_add_a_inv, cs, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}


QDD
qdd_cmult(QDD qdd, uint64_t a, uint64_t N)
{
    // this implements the _controlled_ cmult operation
    // 1. QFT on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);

    // 2. loop over k = {0, n-1}
    uint64_t t = a;
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (BDDVAR i = shor_wires.ctrl_last; i >= shor_wires.ctrl_first; i--) {
        // 2a. double controlled phi_add_mod(a* 2^k)
        cs[1] = i;
        qdd = qdd_phi_add_mod(qdd, cs, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_cmult_inv(QDD qdd, uint64_t a, uint64_t N)
{
    // not quite inverse of above
    // 1. QFT on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT, shor_wires.targ_first, shor_wires.targ_last);
    
    // 2. same loop over k but with phi_add_mod_inv
    uint64_t t = a;
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (BDDVAR i = shor_wires.ctrl_last; i >= shor_wires.ctrl_first; i--) {
        // 2a. double controlled phi_add_mod_inv(a* 2^k)
        cs[1] = i;
        qdd = qdd_phi_add_mod_inv(qdd, cs, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit(qdd, CIRCID_QFT_inv, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_shor_ua(QDD qdd,  uint64_t a, uint64_t N)
{
    LACE_ME;

    // 1. controlled Cmult(a)
    qdd = qdd_cmult(qdd, a, N);

    // 2. controlled swap top/bottom registers
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (uint32_t i = shor_wires.ctrl_first; i <= shor_wires.ctrl_last; i++) {
        qdd = qdd_ccircuit(qdd, CIRCID_swap, cs, i, shor_wires.targ_first+i);
    }

    // 3. controlled Cmult_inv(a^-1)
    uint64_t a_inv = inverse_mod(a, N);
    qdd = qdd_cmult_inv(qdd, a_inv, N);

    return qdd;
}

uint64_t
shor_period_finding(uint64_t a, uint64_t N)
{
    // Circuit (quantum period finding of f(x) = a^x mod N)
    // create QDD |0>|0..001>|0>|0..00>
    uint32_t num_qubits = 2*shor_n + 3;
    bool x[num_qubits];
    for (BDDVAR k = 0; k < num_qubits; k++) x[k] = 0;
    x[shor_wires.ctrl_last] = 1; // set the input reg. to |0...001> = |1>

    QDD qdd = qdd_create_basis_state(num_qubits, x);

    uint64_t as[2*shor_n];
    as[2*shor_n-1] = a;
    uint64_t new_a = a;
    for(int i = 2*shor_n-2; i >= 0; i--) {
        new_a = new_a * new_a;
        new_a = new_a % N;
        as[i] = new_a;
    }

    LACE_ME;

    int m_outcomes[2*shor_n];
    int m_outcome;
    double m_prob;

    for (uint32_t i = 0; i < 2*shor_n; i++) {
        if (testing_mode) assert(qdd_is_unitvector(qdd, num_qubits));

        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);

        // controlled Ua^...
        qdd = qdd_shor_ua(qdd, as[i], N);

        // phase gates based on previous measurement
        int k = 2; // First gate needs to be R^dag(2) = S^dag
        for (int j = i-1; j >= 0; j--) {
            if (m_outcomes[j] == 1)
                qdd = qdd_gate(qdd, GATEID_Rk_dag(k), shor_wires.top);
            k += 1; // R(k) is a (2*pi / 2^k) rotation
        }

        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);

        // measure q0
        //sylvan_clear_cache();
        qdd = qdd_measure_qubit(qdd, shor_wires.top, num_qubits, &m_outcome, &m_prob);
        m_outcomes[i] = m_outcome;

        // make sure q0 is in the |0> state
        if (m_outcome == 1) qdd = qdd_gate(qdd, GATEID_X, shor_wires.top);
        qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
        assert(qddnode_getptrhigh(node) == QDD_TERMINAL);
    }

    // turn measurement outcomes into an integer
    uint64_t res = 0;
    for (uint32_t i = 0; i < 2*shor_n; i++) {
        int index = 2*shor_n-1-i;
        res = (res << 1) + m_outcomes[index];
    }
    return res;
}

void
shor_set_globals(uint64_t a, uint64_t N) 
{
    shor_n = ceil(log2(N)); // number of bits for N (not the number of qubits!)  
    uint64_t p2 = 1;
    for (uint32_t i = 0; i < shor_n; i++) { // LSB in bits[0], MSB in bits[63]
        shor_bits_a[i] = a & p2;
        shor_bits_N[i] = N & p2;
        p2 = p2 << 1;
    }
    // Set wire numbers
    shor_wires.top        = 0;
    shor_wires.ctrl_first = 1;
    shor_wires.ctrl_last  = shor_n;
    shor_wires.helper     = shor_n + 1; // easier to have this in the middle
    shor_wires.targ_first = shor_n + 2;
    shor_wires.targ_last  = 2*shor_n + 2;
}

uint64_t 
my_gcd (uint64_t a, uint64_t b) // clash with gcd in sylvan_mtbdd.c ...
{
  uint64_t c;
  while ( a != 0 ) { c = a; a = b%a;  b = c; }
  return b;
}

uint64_t modpow(uint64_t base, uint64_t exp, uint64_t modulus) {
  base %= modulus;
  uint64_t result = 1;
  while (exp > 0) {
    if (exp & 1) result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp >>= 1;
  }
  return result;
}

uint64_t
shor_post_process(uint64_t N, uint64_t a, uint64_t b, uint64_t denom, bool verbose)
{
    // For b the following is true:
    // b/denom = x/r, where denom = 2^num bits, and r = the period we want.
    // This function tries to find that r.
    // Implementation from [zulehner2018advanced]
    if (b == 0) {
        if (verbose)
            printf("Factorization failed (measured 0)\n");
        return 0;
    }

    int cf_max_size = 100;
    int cf_entries = 0;
    uint64_t cf[cf_max_size];
    uint64_t old_b = b;
	uint64_t old_denom = denom;
	while(b != 0) {
        if (cf_entries >= cf_max_size) {
            printf("please hardcode cf_max_size to something bigger\n"); // I'm sorry
            exit(1);
        }
        cf[cf_entries] = (denom/b);
        cf_entries++;
		uint64_t tmp = denom % b;
		denom = b;
		b = tmp;
	}

    if (verbose) {
        printf("Continued fraction expansion of %ld/%ld = ", b, denom);
        for(int i = 0; i < cf_entries; i++) printf("%ld ", cf[i]);
        printf("\n");
    }

	for(int i=0; i < cf_entries; i++) {
		//determine candidate
		uint64_t denominator = cf[i];
		uint64_t numerator = 1;

		for(int j=i-1; j >= 0; j--) {
			uint64_t tmp = numerator + cf[j]*denominator;
			numerator = denominator;
			denominator = tmp;
		}
        if (verbose)
            printf(" Candidate %ld/%ld: ", numerator, denominator);
		if(denominator > N) {
            if (verbose) {
                printf(" denominator too large (greater than %ld)\n", N);
                printf("Factorization failed\n");
            }
            return 0;
		} else {
			double delta = (double)old_b / (double)old_denom - (double)numerator / (double) denominator;
			if(fabs(delta) < 1.0/(2.0*old_denom)) {
				if(modpow(a, denominator, N) == 1) {
                    if (verbose)
                        printf("found period = %ld\n", denominator);
					if(denominator & 1) {
                        if (verbose)
                            printf("Factorization failed (period is odd)\n");
                        return 0;
					} else {
						uint64_t f1, f2;
						f1 = modpow(a, denominator>>1, N);
						f2 = (f1+1)%N;
						f1 = (f1 == 0) ? N-1 : f1-1;
						f1 = my_gcd(f1, N);
						f2 = my_gcd(f2, N);
                        if (f1 == 1 || f1 == N) {
                            if (verbose)
                                printf("Factorization found trivial factor\n");
                            return 0;
                        }
                        if (verbose) {
                            printf("Factorization succeeded! Non-trivial factors are:\n");
                            printf(" -- gcd(%ld^(%ld/2)-1,%ld)=%ld\n", N, denominator, N, f1);
                            printf(" -- gcd(%ld^(%ld/2)+1,%ld)=%ld\n", N, denominator, N, f2);
                        }
                        return f1;
					}
					break;
				} else {
                    if (verbose)
                        printf("failed\n");
				}
			} else {
                if (verbose)
                    printf("delta is too big (%lf)\n", delta);
			}
		}
	}
    return 0;
}

uint64_t
shor_generate_a(uint64_t N)
{
    uint64_t a;
    do {
        a = rand() % N;
    } while (my_gcd(a, N) != 1 || a == 1);
    return a;
}

uint64_t
run_shor(uint64_t N, uint64_t a, bool verbose)
{
    // The classical part
    if (a == 0) a = shor_generate_a(N);

    shor_set_globals(a, N);
    
    if (verbose) {
        printf("input N        = %ld [", N);
        for (uint32_t i=0; i<shor_n; i++) printf("%d", shor_bits_N[i]);
        printf("]\n");
        printf("n (bits for N) = %d\n",  shor_n);
        printf("random a       = %ld [", a);
        for (uint32_t i=0; i<shor_n; i++) printf("%d", shor_bits_a[i]);
        printf("]\n\n");

        printf("wires:\n");
        printf("top:        %d\n", shor_wires.top);
        printf("ctrl_first: %d\n", shor_wires.ctrl_first);
        printf("ctrl_last:  %d\n", shor_wires.ctrl_last);
        printf("helper:     %d\n", shor_wires.helper);
        printf("targ_first: %d\n", shor_wires.targ_first);
        printf("targ_last:  %d\n\n", shor_wires.targ_last);
    }

    uint64_t b = shor_period_finding(a, N);
    uint64_t denom = 1 << (2*shor_n);

    return shor_post_process(N, a, b, denom, verbose);
}

/******************************</Shor components>******************************/



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

    prob_low  = qdd_unnormed_prob(low,  1, nvars);
    prob_high = qdd_unnormed_prob(high, 1, nvars);
    prob_root = comp_to_prob(comp_value(QDD_AMP(qdd)));
    prob_low  *= prob_root;
    prob_high *= prob_root;
    if (fabs(prob_low + prob_high - 1.0) > TOLERANCE) {
        printf("prob sum = %.5lf (%.5lf + %.5lf)\n", prob_low + prob_high, prob_low, prob_high);
        assert("probabilities don't sum to 1" && false);
    }

    // flip a coin
    float rnd = ((float)rand())/RAND_MAX;
    *m = (rnd < prob_low) ? 0 : 1;
    *p = prob_low;

    // produce post-measurement state
    complex_t norm;
    if (*m == 0) {
        high = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        norm = comp_make(sqrt(prob_low), 0.0);
    }
    else {
        low  = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ZERO);
        norm = comp_make(sqrt(prob_high), 0.0);
    }

    QDD res = qdd_makenode(0, low, high);

    complex_t c = comp_mul(comp_value(QDD_AMP(qdd)), comp_value(QDD_AMP(res)));
    c = comp_div(c, norm);
    AMP new_root_amp = qdd_comp_lookup(c);

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
        prob_roots *= comp_to_prob(comp_value(QDD_AMP(qdd)));
        prob_high = prob_high * prob_roots / prob_path;
        prob_low  = prob_low  * prob_roots / prob_path;

        if (fabs(prob_low + prob_high - 1.0) > TOLERANCE) {
            printf("prob sum = %.10lf\n", prob_low + prob_high);
            assert("probabilities don't sum to 1" && false);
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
        return comp_to_prob(comp_value(QDD_AMP(qdd)));
    }

    // Look in cache
    bool cachenow = 1;
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
    prob_root = comp_to_prob(comp_value(QDD_AMP(qdd)));
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

    complex_t c = comp_one();
    QDD low, high;
    for (;;) {
        c = comp_mul(c, comp_value(QDD_AMP(q)));
        
        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        BDDVAR var     = qddnode_getvar(node);
        qddnode_getchilderen(node, &low, &high);

        // Condition low/high choice on basis state vector[var]
        q = (basis_state[var] == 0) ? low : high;
    }

    return qdd_comp_lookup(c);
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
qdd_stack_matrix(QDD below, BDDVAR k, uint32_t gateid)
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
qdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, uint32_t gateid)
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
qdd_create_single_qubit_gates(BDDVAR n, uint32_t *gateids)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qdd_stack_matrix(prev, k, gateids[k]);
    }
    return prev;
}

QDD
qdd_create_single_qubit_gates_same(BDDVAR n, uint32_t gateid)
{
    // Start at terminal and build backwards
    QDD prev = qdd_bundle_ptr_amp(QDD_TERMINAL, C_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qdd_stack_matrix(prev, k, gateid);
    }
    return prev;
}

QDD
qdd_create_controlled_gate(BDDVAR n, BDDVAR c, BDDVAR t, uint32_t gateid)
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
        ccphase = qdd_bundle_ptr_amp(QDD_PTR(ccphase), C_MIN_ONE);
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
    complex_t c = comp_value(QDD_AMP(qdd));
    c.r = sqrt(c.r*c.r + c.i*c.i);
    c.i = 0.0;
    AMP new_root_amp = qdd_comp_lookup(c);
    QDD res = qdd_bundle_ptr_amp(QDD_PTR(qdd), new_root_amp);
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
            if(!comp_exact_equal(comp_value(amp_a), comp_value(amp_b))){
                if(verbose){
                    _print_bitstring(x, n, true);
                    printf(", amp a ="); comp_print(comp_value(amp_a));
                    printf(" != amp b ="); comp_print(comp_value(amp_b));
                    printf("\n");
                }
                return false;
            }
        }
        else{
            if(!comp_approx_equal(comp_value(amp_a), comp_value(amp_b))){
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
        sum_abs_squares += comp_to_prob(comp_value(a));
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
    return qdd_is_close_to_unitvector(qdd, n, TOLERANCE*10);
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

bool
bit_from_int(uint64_t a, uint8_t index)
{
    // assumes index=0 is the LSB
    uint64_t mask = 1<<index;
    uint64_t res = a & mask;
    res = res>>index;
    return (bool) res;
}

/*******************************<logging stats>********************************/

bool qdd_stats_logging = false;
uint64_t statslog_buffer = 1024;
uint64_t in_buffer = 0;
uint64_t *nodelog;
uint64_t *amp_log;
FILE *qdd_logfile;
uint64_t nodes_peak = 0;
uint64_t logcounter = 0;

void
qdd_stats_start(FILE *out)
{
    if (out == NULL) return;
    qdd_stats_logging = true;
    qdd_logfile = out;
    fprintf(qdd_logfile, "nodes, amps\n");
    nodelog = (uint64_t*) malloc(statslog_buffer * sizeof(uint64_t));
    amp_log = (uint64_t*) malloc(statslog_buffer * sizeof(uint64_t));
    in_buffer = 0;
    nodes_peak = 0;
    logcounter = 0;
}

void
qdd_stats_flush_buffer()
{
    for (uint64_t i = 0; i < in_buffer; i++) {
        fprintf(qdd_logfile, "%ld,%ld\n", nodelog[i], amp_log[i]);
    }
    in_buffer = 0;
}

void
qdd_stats_log(QDD qdd)
{
    if (!qdd_stats_logging) return;

    if (in_buffer >= statslog_buffer) {
        qdd_stats_flush_buffer();
    }

    // Insert info
    uint64_t num_nodes = qdd_countnodes(qdd);
    uint64_t num_amps  = count_amplitude_table_enries();
    nodelog[in_buffer] = num_nodes;
    amp_log[in_buffer] = num_amps;
    in_buffer++;
    logcounter++;

    if (num_nodes > nodes_peak)
        nodes_peak = num_nodes;
}

uint64_t
qdd_stats_get_nodes_peak()
{
    return nodes_peak;
}

uint64_t
qdd_stats_get_logcounter()
{
    return logcounter;
}

void
qdd_stats_finish()
{
    if (!qdd_stats_logging) return;
    qdd_stats_flush_buffer();
    qdd_stats_logging = false;
    nodes_peak = 0;
    logcounter = 0;
    free(nodelog);
    free(amp_log);
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
qdd_fprintdot_label(FILE *out, AMP a)
{
    fprintf(out, ", label=\"");
    if (a == C_ONE) {}
    else if (a == C_ZERO) { fprintf(out, "0"); }
    else {
        complex_t val = comp_value(a);
        if (val.r != 0.0) fprintf(out, "%.3lf", val.r);
        if (val.i > 0.0) fprintf(out, "+%.3lfi", val.i);
        else if (val.i < 0.0) fprintf(out, "%.3lfi", val.i);
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
        qdd_fprintdot_label(out, QDD_AMP(low));
        fprintf(out, "];\n");
    }
    if (draw_zeros || QDD_AMP(high) != C_ZERO) {
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                    QDD_PTR(qdd), QDD_PTR(high));
        qdd_fprintdot_label(out, QDD_AMP(high));
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
    qdd_fprintdot_label(out, QDD_AMP(qdd));
    fprintf(out, "];\n");

    // terminal node
    fprintf(out, "%lu [shape=box, label=\"T\"];\n", QDD_TERMINAL);

    // recursively add nodes
    qdd_fprintdot_rec(out, qdd, draw_zeros);
    qdd_unmark_rec(qdd);

    fprintf(out, "}\n");
}

/**************************</printing & file writing>**************************/
