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


/**
 * QDD node structure (128 bits)
 * 
 * 64 bits low:
 *       4 bits: lower 4 bits of 8 bit variable/qubit number of this node
 *      20 bits: index of edge weight of low edge in Ctable (AMP)
 *      40 bits: low edge pointer to next node (PTR)
 * 64 bits high:
 *       4 bits: upper 4 bits of 8 bit variable/qubit number of this node
 *      20 bits: index of edge weight of high edge in Ctable (AMP)
 *      40 bits: high edge pointer to next node (PTR)
 */
typedef struct __attribute__((packed)) qddnode {
    QDD low, high; // TODO: rename to low, high
} *qddnode_t; // 16 bytes


/**
 * Gets only the AMP information of a QDD edge `q`.
 */
static inline AMP
QDD_AMP(QDD q)
{
    return (q & 0x0fffff0000000000) >> 40; // 20 bits
}

/**
 * Gets only the PTR information of a QDD edge `q`.
 */
static inline PTR
QDD_PTR(QDD q)
{
    return q & 0x000000ffffffffff; // 40 bits
}


/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
qdd_getvar(qddnode_t n)
{
    return (BDDVAR) ( (n->low >> 60) | ((n->high >> 56) & 0xf0) ); // 8 bits
}

/**
 * Gets the low edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qdd_getlow(qddnode_t n)
{
    return (QDD) n->low & 0x0fffffffffffffff; // 60 bits
}

/**
 * Gets the high edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qdd_gethigh(qddnode_t n)
{
    return (QDD) n->high & 0x0fffffffffffffff; // 60 bits
}

/**
 * Gets only the PTR of the low edge of `n`.
 */
static inline PTR
qdd_getptrlow(qddnode_t n)
{
    return (PTR) QDD_PTR(n->low);
}

/**
 * Gets only the PTR of the high edge of `n`.
 */
static inline PTR
qdd_getptrhigh(qddnode_t n)
{
    return (PTR) QDD_PTR(n->high);
}

/**
 * Gets only the AMP of the low edge of `n`.
 */
static inline AMP
qdd_getamplow(qddnode_t n)
{
    return (AMP) QDD_AMP(n->low);
}

/**
 * Gets only the AMP of the high edge of `n`.
 */
static inline AMP
qdd_getamphigh(qddnode_t n)
{
    return (AMP) QDD_AMP(n->high);
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

/**
 * Gets the node `p` is pointing to.
 */
static inline qddnode_t
QDD_GETNODE(PTR p)
{
    return (qddnode_t) llmsset_index_to_ptr(nodes, p);
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
    n->low  = qdd_bundle_low(var, low, a);
    n->high = qdd_bundle_high(var, high, b);
}



static PTR
//_qdd_makenode(BDDVAR var, QDD low, QDD high)
_qdd_makenode(BDDVAR var, PTR low, PTR high, AMP a, AMP b)
{
    struct qddnode n;

    qddnode_make(&n, var, low, high, a, b);

    PTR result;
    int created;
    PTR index = llmsset_lookup(nodes, n.low, n.high, &created);
    if (index == 0) {
        LACE_ME;

        mtbdd_refs_push(n.low);//mtbdd_refs_push(low);
        mtbdd_refs_push(n.high);//mtbdd_refs_push(high);
        sylvan_gc();
        mtbdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.low, n.high, &created);
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

static inline QDD // (PTR and AMP, but the amp is the norm weight from below)
qdd_makenode(BDDVAR var, QDD low_edge, QDD high_edge)
{ 
    PTR low_ptr  = QDD_PTR(low_edge);
    AMP low_amp  = QDD_AMP(low_edge);
    PTR high_ptr = QDD_PTR(high_edge);
    AMP high_amp = QDD_AMP(high_edge);

    PTR p;
    AMP norm;

    // Edges with weight 0 should point straight to terminal.
    if(low_amp  == C_ZERO) low_ptr  = QDD_TERMINAL;
    if(high_amp == C_ZERO) high_ptr = QDD_TERMINAL;

    // If both edges are the same, normalize
    if(low_ptr == high_ptr && low_amp == high_amp) {
        // Both 0 (1), norm factor is 0 (1)
        if(low_amp == C_ZERO || low_amp == C_ONE){
            return qdd_bundle_ptr_amp(low_ptr, low_amp);
        }
        else{
            // Both amps are equal, but not 0 or 1
            norm = low_amp;
            p = _qdd_makenode(var, low_ptr, high_ptr, C_ONE, C_ONE);
            return qdd_bundle_ptr_amp(p, norm);
        }
    }
    else{
        // If the edges are not the same
        if(low_amp != C_ZERO){
            // Normalize using low
            norm     = low_amp;
            low_amp  = C_ONE;
            high_amp = Cdiv(high_amp, norm);
        }
        else{
            // Unless low amp is 0, then norm using high
            norm     = high_amp;
            high_amp = C_ONE;
        }
        p = _qdd_makenode(var, low_ptr, high_ptr, low_amp, high_amp);
        return qdd_bundle_ptr_amp(p, norm);
    }
}


// TODO: also return and deal with normalization parameter
QDD qdd_apply_gate(QDD q, uint32_t gate, BDDVAR qubit)
{
    if(QDD_PTR(q) == QDD_TERMINAL)
        return q;
    
    qddnode_t node = QDD_GETNODE(QDD_PTR(q));

    BDDVAR var = qdd_getvar(node);

    pprint_qddnode(node);
    
    // "above" the desired qubit in the QDD
    if(var < qubit){
        QDD low  = qdd_apply_gate(node->low,  gate, qubit);
        QDD high = qdd_apply_gate(node->high, gate, qubit);
        return qdd_makenode(var, low, high);
    }
    // exactly at or passed qubit
    else {
        PTR l, h;
        AMP a, b;
        if(var == qubit){
            l = qdd_getptrlow(node);
            h = qdd_getptrhigh(node);
            a = qdd_getamplow(node);
            b = qdd_getamphigh(node);
        }
        else {
            l = QDD_PTR(q);
            h = QDD_PTR(q);
            a = C_ONE;
            b = C_ONE;
        }
        // E.g. X gate:
        AMP a_new = b; // TODO: actual matrix multiplication
        AMP b_new = a; // TODO: actual matrix multiplication
        QDD low  = qdd_bundle_ptr_amp(l, a_new);
        QDD high = qdd_bundle_ptr_amp(h, b_new);
        return qdd_makenode(var, low, high);
    }
}


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


    qddnode_t node_a = QDD_GETNODE(QDD_PTR(a));
    qddnode_t node_b = QDD_GETNODE(QDD_PTR(b));

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
        aLow = qdd_getlow(node_a);
        aHigh = qdd_gethigh(node_b);
    }
    if (level == vb) {
        bLow = qdd_getlow(node_b);
        bHigh = qdd_gethigh(node_b);
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

    //result = qdd_makenode(level, low, high);
    result = QDD_PTR(qdd_makenode(level, low, high));

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
            qddnode_t qn = QDD_GETNODE(QDD_PTR(q));
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
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        BDDVAR var     = qdd_getvar(node);
        pprint_qddnode(node);

        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // Condition low/high choice on basis state vector[var]
        if (basis_state[var] == 0)
            q = node->low;
        else
            q = node->high;
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

        low_child = QDD_PTR(qdd_makenode(k, n.low, n.high));
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


