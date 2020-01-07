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


// I don't know why this is defined in the source file, and AMP in the header,
// however moving it to the header file creates issues with IntelliSense (but 
// not for the actual compiler though).
//typedef uint64_t        PTR; // pointer


static inline AMP __attribute__((unused))
QDD_AMP(QDD q)
{
    // Mask out the 4 top bits (var) and shift right
    return (q & 0x0fffff0000000000) >> 40; // 20 bits
}

static inline PTR __attribute__((unused))
QDD_PTR(QDD q)
{
    return q & 0x000000ffffffffff; // 40 bits
}

/**
 * QDD node structure
 *
 */
typedef struct __attribute__((packed)) qddnode {
    //QDD a, b;
    uint64_t low : 40;
    uint64_t high : 40;
    uint64_t amp_low : 20;
    uint64_t amp_high : 20;
    uint64_t var : 8;
} *qddnode_t; // 16 bytes

static inline void pprint_qddnode(qddnode_t n)
{
    printf("[var=%d, low=%p, high=%p, a=%p, b=%p]\n", 
             n->var, n->low, n->high, n->amp_low, n->amp_high);
}

static inline QDD bundle_edge(PTR p, AMP a)
{
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= 0x00000000000fffff);
    return a << 40 | p;
}

static inline qddnode_t
QDD_GETNODE(QDD q)
{
    return (qddnode_t) llmsset_index_to_ptr(nodes, QDD_PTR(q));
}

//static inline QDD __attribute__((unused))
//qdd_getlow(qddnode_t n)
//{
//    return n->a & 0x0fffffffffffffff;
//}
static inline QDD __attribute__((unused))
qdd_getlow_bundled(qddnode_t n)
{
    return bundle_edge(n->low, n->amp_low);
}

//static inline QDD __attribute__((unused))
//qdd_gethigh(qddnode_t n)
//{
//    return n->b & 0x0fffffffffffffff;
//}
static inline QDD __attribute__((unused))
qdd_gethigh_bundled(qddnode_t n)
{
    return bundle_edge(n->high, n->amp_high);
}

//static inline BDDVAR __attribute__((unused))
//qdd_getvar(qddnode_t n)
//{
//    return (BDDVAR) ( (n->a >> 60) | ((n->b >> 56) & 0xf0) );
//}

//static inline QDD qdd_make(PTR p, AMP a)
//{
//    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
//    assert (a <= 0x00000000000fffff);
//    return a << 40 | p;
//} // replaced with bundle_edge(PTR p, AMP a)

//static inline void qddnode_make(qddnode_t n, BDDVAR var, QDD low, QDD high)
//{
//    assert (var <= 0xff);
//    // changed n->a and n->b arround (maybe better to call qdd_getlow/high)
//    n->a = low   |  (((uint64_t)var) << 60);
//    n->b = high  | ((((uint64_t)var) << 56) & 0xf000000000000000);
//}
static inline void 
qddnode_make(qddnode_t n, BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    assert (var <= 0xff);
    n->low      = low;
    n->high     = high;
    n->amp_low  = a;
    n->amp_high = b;
    n->var      = var;
}


/*
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
} */
static PTR
_qdd_makenode(BDDVAR var, PTR low_node, PTR high_node, AMP a, AMP b)
{
    struct qddnode n;

    qddnode_make(&n, var, low_node, high_node, a, b);

    QDD low_edge;
    QDD high_edge;

    PTR result;
    int created;
    
    // For the llmsset_lookup function, the pointers to low and high also need
    // to contain the information about the amplitudes as well for this function
    // to work correctly.
    low_edge  = bundle_edge(n.low, n.amp_low);
    high_edge = bundle_edge(n.high, n.amp_high);
    PTR index = llmsset_lookup(nodes, low_edge, high_edge, &created);
    if (index == 0) {
        LACE_ME;

        low_edge  = bundle_edge(low_node, a);
        high_edge = bundle_edge(low_node, b);
        mtbdd_refs_push(low_edge);
        mtbdd_refs_push(high_edge);
        sylvan_gc();
        mtbdd_refs_pop(2);

        low_edge  = bundle_edge(n.low, n.amp_low);
        high_edge = bundle_edge(n.high, n.amp_high);
        index = llmsset_lookup(nodes, low_edge, high_edge, &created);
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

//static inline PTR qdd_makenode(BDDVAR var, QDD low, QDD high)
//{
//    
//    return low == high ? low : _qdd_makenode(var, low, high);
//}
static inline PTR 
qdd_makenode(BDDVAR var, QDD low_edge, QDD high_edge)
{
    if(low_edge == high_edge)
        return QDD_PTR(low_edge);
    else
        return _qdd_makenode(var, QDD_PTR(low_edge), QDD_PTR(high_edge), 
                                  QDD_AMP(low_edge), QDD_AMP(high_edge));
}

//static inline QDD qdd_makenode2(BDDVAR var, AMP a, QDD low, QDD high)
//{
//    return qdd_make(qdd_makenode(var, low,high), a);
//}

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

    BDDVAR va = node_a->var; //BDDVAR va = qdd_getvar(na);
    BDDVAR vb = node_b->var; //BDDVAR vb = qdd_getvar(nb);
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
        aLow = qdd_getlow_bundled(node_a);
        aHigh = qdd_gethigh_bundled(node_a);
    }
    if (level == vb) {
        bLow = qdd_getlow_bundled(node_b);
        bHigh = qdd_gethigh_bundled(node_b);
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
            if (qn->var == mtbddnode_getvariable(n_vars)) {

                q = *str ? qdd_gethigh_bundled(qn) : qdd_getlow_bundled(qn);
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
    // TODO: actually store amps in a table and use AMP as an index to that
    //       table (also for create_all_zero_state())
    // TODO: clean up this function
    if (q == sylvan_false) return 0;

    AMP a = 1;
    for (;;) {
        // multiply `a` with amplitude of current QDD edge
        //a = Cmul(QDD_AMP(q), a);

        // temp:
        printf("amp=%d ", QDD_AMP(q));
        a = a * QDD_AMP(q);
        
        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(q);
        BDDVAR var     = node->var;
        printf("at node with var=%d\n", var);

        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_PTR(QDD_TERMINAL)) break;

        // Condition low/high on bosis state vector[var]
        if (basis_state[var] == 0)
            q = qdd_getlow_bundled(node);
        else
            q = qdd_gethigh_bundled(node);
    }

    return a;
}

QDD
create_all_zero_state(int n_qubits)
{
    assert(n_qubits >= 1);

    struct qddnode n;
    
    // start at terminal, and build backwards
    PTR prev_node = QDD_TERMINAL;
    for(int k = n_qubits-1; k >= 0; k--){
        
        qddnode_make(&n, k, prev_node, QDD_TERMINAL, 1, NIL);
        // TODO: replace amps with call to amp table

        printf("Created the following node:\n");
        pprint_qddnode(&n);

        prev_node = qdd_makenode(k, qdd_getlow_bundled(&n), qdd_gethigh_bundled(&n));
        printf("With index in the nodetable = %p\n", prev_node);
    }

    QDD root_edge = bundle_edge(prev_node, 1);
    return root_edge;
}

// just for testing TODO: remove
void
init_amplitude_table()
{
    qdd_complex_init();
}


