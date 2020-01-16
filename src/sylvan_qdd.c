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

#include <sylvan_qdd_int.h>

static int granularity = 1; // default




static inline AMP
QDD_AMP(QDD q)
{
    // Mask out the 4 top bits (var) and shift right
    return (q & 0x0fffff0000000000) >> 40; // 20 bits
}

static inline PTR
QDD_PTR(QDD q)
{
    return q & 0x000000ffffffffff; // 40 bits
}

/**
 * QDD node structure
 * 
 * TODO: document 
 *
 */
typedef struct __attribute__((packed)) qddnode {
    QDD a, b; // TODO: rename to low, high
} *qddnode_t; // 16 bytes


static inline BDDVAR
qdd_getvar(qddnode_t n)
{
    // 4 highest bits of `b` followed by the 4 highest bits of `a`.
    return (BDDVAR) ( (n->a >> 60) | ((n->b >> 56) & 0xf0) ); // 8 bits
}

/**
 * Gets the low edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qdd_getlow(qddnode_t n)
{
    return (QDD) n->a & 0x0fffffffffffffff;
}

/**
 * Gets the high edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qdd_gethigh(qddnode_t n)
{
    return (QDD) n->b & 0x0fffffffffffffff;
}

/**
 * Gets only the PTR of the low edge of `n`.
 */
static inline PTR
qdd_getptrlow(qddnode_t n)
{
    return (PTR) QDD_PTR(n->a);
}

/**
 * Gets only the PTR of the high edge of `n`.
 */
static inline PTR
qdd_getptrhigh(qddnode_t n)
{
    return (PTR) QDD_PTR(n->b);
}

/**
 * Gets only the AMP of the low edge of `n`.
 */
static inline AMP
qdd_getamplow(qddnode_t n)
{
    return (AMP) QDD_AMP(n->a);
}

/**
 * Gets only the AMP of the high edge of `n`.
 */
static inline AMP
qdd_getamphigh(qddnode_t n)
{
    return (AMP) QDD_AMP(n->b);
}

/**
 * Pretty prints the information contained in `n`.
 */
static inline void pprint_qddnode(qddnode_t n)
{
    printf("[var=%d, low=%p, high=%p, a=%p, b=%p]\n", 
             qdd_getvar(n),
             qdd_getptrlow(n),
             qdd_getptrhigh(n),
             qdd_getamplow(n),
             qdd_getamphigh(n));
}

static inline qddnode_t
QDD_GETNODE(QDD q)
{
    // TODO: I think it's better if this function takes a PTR as argument,
    // rather than hiding the QDD_PTR(q) from the outside.
    return (qddnode_t) llmsset_index_to_ptr(nodes, QDD_PTR(q));
}


static inline QDD
qdd_bundle_ptr_amp(PTR p, AMP a)
{
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= 0x00000000000fffff);
    return (a << 40 | p);
}


static inline QDD
qdd_bundle_low(BDDVAR var, PTR p, AMP a)
{
    // on the low edge we store the bottom 4 bits of the 8 bit var
    assert (var <= 0xff);
    QDD q = qdd_bundle_ptr_amp(p, a);
    q = ((uint64_t)var << 60) | q;
    return q;
}


static inline QDD
qdd_bundle_high(BDDVAR var, PTR p, AMP a)
{
    // on the high edge we store the top 4 bits of the 8 bit var
    assert (var <= 0xff);
    QDD q = qdd_bundle_ptr_amp(p, a);
    return (((uint64_t)var & 0xf0) << 56) | q;
}


static inline void 
qddnode_make(qddnode_t n, BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    n->a = qdd_bundle_low(var, low, a);
    n->b = qdd_bundle_high(var, high, b);
}



static PTR
//_qdd_makenode(BDDVAR var, QDD low, QDD high)
_qdd_makenode(BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    struct qddnode n;

    qddnode_make(&n, var, low, high, a, b);

    PTR result;
    int created;
    PTR index = llmsset_lookup(nodes, n.a, n.b, &created);
    if (index == 0) {
        LACE_ME;

        mtbdd_refs_push(n.a);//mtbdd_refs_push(low);
        mtbdd_refs_push(n.b);//mtbdd_refs_push(high);
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

static inline PTR 
qdd_makenode(BDDVAR var, QDD low_edge, QDD high_edge)
{
    if(low_edge == high_edge)
        return QDD_PTR(low_edge);
    else
        return _qdd_makenode(var, QDD_PTR(low_edge), QDD_PTR(high_edge), 
                                  QDD_AMP(low_edge), QDD_AMP(high_edge));
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


    qddnode_t node_a = QDD_GETNODE(a);
    qddnode_t node_b = QDD_GETNODE(b);

    BDDVAR va = qdd_getvar(node_a);
    BDDVAR vb = qdd_getvar(node_b);
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
        // TODO: fix
        //aLow = qdd_getlow_bundled(node_a);
        //aHigh = qdd_gethigh_bundled(node_a);
    }
    if (level == vb) {
        // TODO: fix
        //bLow = qdd_getlow_bundled(node_b);
        //bHigh = qdd_gethigh_bundled(node_b);
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

        a = Cmul(QDD_AMP(q), a);
        *str = (rand() & 0x2000) == 0;

        if (q != QDD_ONE) {
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

AMP
qdd_get_amplitude(QDD q, bool* basis_state)
{
    if (q == sylvan_false) return 0;

    AMP a = C_ONE;
    for (;;) {
        a = Cmul(a, QDD_AMP(q));

        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(q);
        BDDVAR var     = qdd_getvar(node);
        pprint_qddnode(node);

        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // Condition low/high choice on basis state vector[var]
        if (basis_state[var] == 0)
            q = node->a;
        else
            q = node->b;
    }

    printf("amplitude in Ctable:");
    Cprint(Cvalue(a));
    printf("\n");
    // TODO: return complex struct instead of the index.
    return a;
}

QDD
create_all_zero_state(int n_qubits)
{
    assert(n_qubits >= 1);

    struct qddnode n;
    
    // start at terminal, and build backwards
    PTR low_child = QDD_TERMINAL;
    for(int k = n_qubits-1; k >= 0; k--){
        
        qddnode_make(&n, k, low_child, QDD_TERMINAL, C_ONE, C_ZERO);
        // is it correct to use CONE and CZRO here?

        printf("Created the following node:\n");
        pprint_qddnode(&n);

        low_child = qdd_makenode(k, n.a, n.b);
        pprint_qddnode(QDD_GETNODE(low_child));
        printf("With index in the nodetable = %p\n", low_child);
    }

    QDD root_edge = qdd_bundle_ptr_amp(low_child, C_ONE);
    return root_edge;
}

// just for testing TODO: do this somewhere better
void
init_amplitude_table()
{
    qdd_complex_init();
}


