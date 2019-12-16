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

#include <sylvan_int.h>
#include <sylvan_qdd.h>

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


static int granularity = 1; // default


typedef uint64_t        PTR; // pointer


static inline AMP __attribute__((unused))
QDD_AMP(QDD q)
{
    return q >> 40; // 40 bits
}

static inline PTR __attribute__((unused))
QDD_PTR(QDD q)
{
    return q & 0x000000ffffffffff; // 40 bits
}

/**
 * BDD/MTBDD node structure
 *
 *
 */
typedef struct __attribute__((packed)) qddnode {
    QDD a, b;
} *qddnode_t; // 16 bytes

static inline qddnode_t
QDD_GETNODE(QDD q)
{
    return (qddnode_t) llmsset_index_to_ptr(nodes, QDD_PTR(q));
}

static inline QDD __attribute__((unused))
qdd_getlow(qddnode_t n)
{
    return n->b & 0x0fffffffffffffff;
}

static inline QDD __attribute__((unused))
qdd_gethigh(qddnode_t n)
{
    return n->b & 0x0fffffffffffffff;
}

static inline BDDVAR __attribute__((unused))
qdd_getvar(qddnode_t n)
{
    return (BDDVAR) ( (n->b >> 60) | ((n->b >> 56) & 0xf0) );
}

static inline QDD qdd_make(PTR p, AMP a)
{
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= 0x00000000000fffff);
    return a << 40 | p;
}

static inline void qddnode_make(qddnode_t n, BDDVAR var, QDD low, QDD high)
{
    assert (var <= 0xff);
    n->b = low   |  (((uint64_t)var) << 60);
    n->a = high  | ((((uint64_t)var) << 56) & 0xf000000000000000);
}


static PTR
_qdd_makenode(BDDVAR var, QDD low, QDD high)
{
    struct qddnode n;

    qddnode_make(&n, var, low, high);

    PTR result;
    int created;
    PTR index = llmsset_lookup(nodes, n.a, n.b, &created);
    if (index == 0) {
        LACE_ME;

        mtbdd_refs_push(low);
        mtbdd_refs_push(high);
        sylvan_gc();
        mtbdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.a, n.b, &created);
        if (index == 0) {
            fprintf(stderr, "BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(BDD_NODES_CREATED);
    else sylvan_stats_count(BDD_NODES_REUSED);

    result = index;
    return result;
}

static inline PTR qdd_makenode(BDDVAR var, QDD low, QDD high)
{
    return low == high ? low : _qdd_makenode(var, low, high);
}

static inline QDD qdd_makenode2(BDDVAR var, AMP a, QDD low, QDD high)
{
    return qdd_make(qdd_makenode(var, low,high), a);
}

//TODO: implement a normalizing make_node code


/**
 * Implementation of plus
 */
TASK_IMPL_3(BDD, qdd_plus, QDD, a, QDD, b, BDDVAR, prev_level)
{
    //AMP aa = QDD_AMP(a);
    //AMP ab = QDD_AMP(b);

    /* Terminal cases */
    if (a == sylvan_true) return b;
    if (b == sylvan_true) return a;
    if (a == sylvan_false) return sylvan_false;
    if (b == sylvan_false) return sylvan_false;
    if (a == b) return a;

    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(BDD_AND);


    qddnode_t na = QDD_GETNODE(a);
    qddnode_t nb = QDD_GETNODE(b);

    BDDVAR va = qdd_getvar(na);
    BDDVAR vb = qdd_getvar(nb);
    BDDVAR level = va < vb ? va : vb;

    int cachenow = granularity < 2 || prev_level == 0 ? 1 : prev_level / granularity != level / granularity;
    if (cachenow) {
        BDD result;
        if (cache_get3(CACHE_BDD_AND, a, b, sylvan_false, &result)) {
            sylvan_stats_count(BDD_AND_CACHED);
            return result;
        }
    }

    // Get cofactors
    QDD aLow = a, aHigh = a;
    QDD bLow = b, bHigh = b;
    if (level == va) {
        aLow = qdd_getlow(na);
        aHigh = qdd_gethigh(na);
    }
    if (level == vb) {
        bLow = qdd_getlow(nb);
        bHigh = qdd_gethigh(nb);
    }

    // Recursive computation
    QDD low = sylvan_invalid,
       high = sylvan_invalid,
     result;

    int n=0;

    if (aHigh == sylvan_true) {
        high = bHigh;
    } else if (aHigh == sylvan_false || bHigh == sylvan_false) {
        high = sylvan_false;
    } else if (bHigh == sylvan_true) {
        high = aHigh;
    } else {
        bdd_refs_spawn(SPAWN(qdd_plus, aHigh, bHigh, level));
        n=1;
    }

    if (aLow == sylvan_true) {
        low = bLow;
    } else if (aLow == sylvan_false || bLow == sylvan_false) {
        low = sylvan_false;
    } else if (bLow == sylvan_true) {
        low = aLow;
    } else {
        low = CALL(qdd_plus, aLow, bLow, level);
    }

    if (n) {
        bdd_refs_push(low);
        high = bdd_refs_sync(SYNC(qdd_plus));
        bdd_refs_pop(1);
    }

    result = qdd_makenode(level, low, high);

    if (cachenow) {
        if (cache_put3(CACHE_BDD_AND, a, b, sylvan_false, result)) sylvan_stats_count(BDD_AND_CACHEDPUT);
    }

    return result;
}

AMP
qdd_sample(QDD q, BDDVAR vars, bool* str)
{
    if (q == sylvan_false) return 0;
    if (sylvan_set_isempty(vars)) return 1;

    AMP a = 1;
    for (;;) {
        mtbddnode_t n_vars = MTBDD_GETNODE(vars);

        a *= QDD_AMP(q);
        *str = (rand() & 0x2000) == 0;

        if (q != ONE) {
            qddnode_t qn = QDD_GETNODE(q);
            if (qdd_getvar(qn) == mtbddnode_getvariable(n_vars)) {

                q = *str ? qdd_gethigh(qn) : qdd_getlow(qn);
            }
        }

        vars = node_high(vars, n_vars);
        if (sylvan_set_isempty(vars)) break;
        str++;
    }

    return 1;
}


