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

#include "sylvan_qdd_int.h"

//static int granularity = 1; // default

// TODO: make a cleaner separation between the code in
// - sylvan_qdd.h
// - sylvan_qdd.c
// - sylvan_qdd_int.h
// - sylvan_qdd_int.c
// consistent with for example the equivalent mtbdd files.

/**
 * QDD node structure (128 bits)
 * 
 * 64 bits low:
 *       4 bits: lower 4 bits of 8 bit variable/qubit number of this node
 *       1 bit:  marked/unmarked flag
 *      19 bits: index of edge weight of low edge in Ctable (AMP)
 *      40 bits: low edge pointer to next node (PTR)
 * 64 bits high:
 *       4 bits: upper 4 bits of 8 bit variable/qubit number of this node
 *       1 bit:  (unused)
 *      19 bits: index of edge weight of high edge in Ctable (AMP)
 *      40 bits: high edge pointer to next node (PTR)
 */
typedef struct __attribute__((packed)) qddnode {
    QDD low, high;
} *qddnode_t; // 16 bytes


/**
 * Gets only the AMP information of a QDD edge `q`.
 */
static inline AMP
QDD_AMP(QDD q)
{
    return (q & 0x07ffff0000000000) >> 40; // 19 bits
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
 *  For caching, we need to uniquely identify gates (and which qubit they are 
 *  applied on)
 *  Uses lowest 40 bits (and should not use more):
 *  [      24 bits blank     | 8 bit c | 8 bit t |     24 bits gateid     ]
 *  Set control = 0 for single qubit gates
 *  (maybe put this elsewhere?)
 */
static inline uint64_t
GATE_OPID(uint32_t gateid, BDDVAR control, BDDVAR target)
{
    uint64_t c = control;
    uint64_t t = target;
    uint64_t res = c<<32 | t<<24 | gateid;
    return res;
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
    return (QDD) n->low & 0x07ffffffffffffff; // 59 bits
}

/**
 * Gets the high edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qdd_gethigh(qddnode_t n)
{
    return (QDD) n->high & 0x07ffffffffffffff; // 59 bits
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
 * Gets the value of the "marked" flag.
 */
static inline bool
qddnode_getmark(qddnode_t n)
{
    return n->low & 0x0800000000000000 ? 1 : 0;
}

/**
 * Sets the value of the "marked" flag to `mark`.
 */
static inline void
qddnode_setmark(qddnode_t n, bool mark)
{
    if (mark) n->low |= 0x0800000000000000; // set 5th bit from left to 1
    else n->low &= 0xf7ffffffffffffff;      // set 5th bit from left to 0
}

/**
 * Pretty prints the information contained in `n`.
 */
static void pprint_qddnode(qddnode_t n)
{
    AMP amp_low  = qdd_getamplow(n);
    AMP amp_high = qdd_getamphigh(n);
    printf("[var=%d, low=%lx, high=%lx, ", 
             qdd_getvar(n),
             qdd_getptrlow(n),
             qdd_getptrhigh(n));
    if(amp_low == C_ZERO)      printf("a=C_ZERO,  ");
    else if(amp_low == C_ONE)  printf("a=C_ONE,   ");
    else                       printf("a=%lx, ",amp_low);
    if(amp_high == C_ZERO)     printf("b=C_ZERO ");
    else if(amp_high == C_ONE) printf("b=C_ONE, ");
    else                       printf("b=%lx", amp_high);
    printf("]\n");
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


static inline QDD
qdd_bundle_ptr_amp(PTR p, AMP a)
{
    assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
    assert (a <= 0x000000000007ffff);
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

    // If both low and high are the same (both PTR and AMP) return low
    if(low_ptr == high_ptr && low_amp == high_amp) {
        return low_edge;
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
    qdd_unmark_rec(qdd_getlow(n));
    qdd_unmark_rec(qdd_gethigh(n));
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
    return 1 + qdd_nodecount_mark(qdd_getlow(n)) + qdd_nodecount_mark(qdd_gethigh(n));
}

uint64_t
qdd_countnodes(QDD qdd)
{
    uint64_t res = qdd_nodecount_mark(qdd) + 1; // (+ 1 for terminal "node")
    qdd_unmark_rec(qdd);
    return res;
}


TASK_IMPL_3(QDD, qdd_gate, QDD, q, uint32_t, gate, BDDVAR, qubit)
{
    if(QDD_PTR(q) == QDD_TERMINAL){
        // Should do the same as 'passed qubit' (i.e. consider terminal as 
        // a node with very high var)
        // (TODO: maybe handle this in QDD_GETNODE)
        PTR l, h;
        AMP a, b;
        l = QDD_PTR(q);
        h = QDD_PTR(q);
        a = C_ONE;
        b = C_ONE;

        AMP a_u00 = Cmul(a, gates[gate][0]);
        AMP a_u10 = Cmul(a, gates[gate][2]);
        AMP b_u01 = Cmul(b, gates[gate][1]);
        AMP b_u11 = Cmul(b, gates[gate][3]);

        QDD qdd1 = qdd_makenode(qubit, qdd_bundle_ptr_amp(l, a_u00), 
                                       qdd_bundle_ptr_amp(l, a_u10));
        QDD qdd2 = qdd_makenode(qubit, qdd_bundle_ptr_amp(h, b_u01),
                                       qdd_bundle_ptr_amp(h, b_u11));
        
        QDD sum = qdd_plus(qdd1, qdd2);
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(sum));
        QDD res = qdd_bundle_ptr_amp(QDD_PTR(sum), new_root_amp);
        return res;
    }
    
    // get node info
    qddnode_t node = QDD_GETNODE(QDD_PTR(q));
    BDDVAR var = qdd_getvar(node);
    
    // "above" the desired qubit in the QDD
    if(var < qubit){

        bool cachenow = 1;
        if (cachenow) {
            QDD res;
            // check if this calculation has already been done before for this node/gate
            if (cache_get3(CACHE_QDD_GATE, GATE_OPID(gate, 0, qubit), q, sylvan_false, &res)) {
                sylvan_stats_count(QDD_GATE_CACHED);
                return res;
            }
        }
        
        QDD res_low, res_high;
        bdd_refs_spawn(SPAWN(qdd_gate, node->high, gate, qubit));
        res_low = CALL(qdd_gate, node->low, gate, qubit);
        bdd_refs_push(res_low);
        res_high = bdd_refs_sync(SYNC(qdd_gate));
        bdd_refs_pop(1);

        QDD res  = qdd_makenode(var, res_low, res_high);
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

        if (cachenow) {
            if (cache_put3(CACHE_QDD_GATE, GATE_OPID(gate, 0, qubit), q, sylvan_false, res)) sylvan_stats_count(QDD_GATE_CACHEDPUT);
        }
        return res;
    }
    else {
        PTR l, h;
        AMP a, b;
        // exactly at qubit
        if(var == qubit){
            l = qdd_getptrlow(node);
            h = qdd_getptrhigh(node);
            a = qdd_getamplow(node);
            b = qdd_getamphigh(node);
        }
        // passed qubit (node with var == qubit was a don't-care)
        else {
            l = QDD_PTR(q);
            h = QDD_PTR(q);
            a = C_ONE;
            b = C_ONE;
        }

        AMP a_u00 = Cmul(a, gates[gate][0]);
        AMP a_u10 = Cmul(a, gates[gate][2]);
        AMP b_u01 = Cmul(b, gates[gate][1]);
        AMP b_u11 = Cmul(b, gates[gate][3]);

        QDD qdd1 = qdd_makenode(qubit, qdd_bundle_ptr_amp(l, a_u00), 
                                       qdd_bundle_ptr_amp(l, a_u10));
        QDD qdd2 = qdd_makenode(qubit, qdd_bundle_ptr_amp(h, b_u01),
                                       qdd_bundle_ptr_amp(h, b_u11));
        QDD sum = qdd_plus(qdd1, qdd2);

        // multiply root amp of sum with input root amp
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(sum));
        QDD res = qdd_bundle_ptr_amp(QDD_PTR(sum), new_root_amp);
        return res;
    }
}

TASK_IMPL_4(QDD, qdd_cgate, QDD, q, uint32_t, gate, BDDVAR, c, BDDVAR, t)
{
    // recursively look for control qubit in QDD, call single qubit gate
    // function on high edge of node with var == c
    assert(c < t);

    if(QDD_PTR(q) == QDD_TERMINAL){
        // passed qubit (node with var == c was a don't-care)
        QDD res_low, res_high;
        res_low  = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        res_high = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        res_high = qdd_gate(res_high, gate, t);
        QDD res = qdd_makenode(c, res_low, res_high);
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
        return res;
    }

    // get node info
    qddnode_t node = QDD_GETNODE(QDD_PTR(q));
    BDDVAR var = qdd_getvar(node);

    // "above" the desired qubit in the QDD (this is where the recursive stuff
    // of cgate happens, once the control qubit has been reached, the recursive
    // search for target qubit continues in qdd_gate function)
    if(var < c){
        QDD res_low, res_high;

        bool cachenow = 1;
        if (cachenow) {
            QDD res;
            // check if this calculation has already been done before for this node/gate
            if (cache_get3(CACHE_QDD_CGATE, GATE_OPID(gate, c, t), q, sylvan_false, &res)) {
                sylvan_stats_count(QDD_CGATE_CACHED);
                return res;
            }
        }

        bdd_refs_spawn(SPAWN(qdd_cgate, node->high, gate, c, t));
        res_low = CALL(qdd_cgate, node->low, gate, c, t);
        bdd_refs_push(res_low);
        res_high = bdd_refs_sync(SYNC(qdd_gate));
        bdd_refs_pop(1);

        QDD res  = qdd_makenode(var, res_low, res_high);
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

        if (cachenow) {
            if (cache_put3(CACHE_QDD_CGATE, GATE_OPID(gate, c, t), q, sylvan_false, res)) sylvan_stats_count(QDD_GATE_CACHEDPUT);
        }
        return res;
    }
    else {
        QDD res_low, res_high;
        //PTR l, h;
        //AMP a, b;
        // exactly at qubit 'c'
        if(var == c){
            res_low = node->low;
            res_high = node->high;
        }
        // passed qubit (node with var == c was a don't-care)
        else {
            res_low  = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
            res_high = qdd_bundle_ptr_amp(QDD_PTR(q), C_ONE);
        }

        // apply gate to branch where 'c' == 1, leave low branch unchanged
        res_high = qdd_gate(res_high, gate, t);
        QDD res = qdd_makenode(c, res_low, res_high); 
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
        return res;
    }
}

TASK_IMPL_2(QDD, qdd_plus, QDD, a, QDD, b)
{
    AMP amp_a = QDD_AMP(a);
    AMP amp_b = QDD_AMP(b);
    
    // Optimization base/terminal cases
    // (are not required for this function to function correctly,
    // but can save work)
    if(amp_a == C_ZERO) return b;
    if(amp_b == C_ZERO) return a;

    // Get info from node (unless terminal)
    BDDVAR var_a, var_b;
    QDD a_low, a_high;
    QDD b_low, b_high;
    if(QDD_PTR(a) == QDD_TERMINAL){
        var_a = 9000; // TODO: do this slightly cleaner
    }
    else{
        qddnode_t node_a = QDD_GETNODE(QDD_PTR(a));
        var_a  = qdd_getvar(node_a);
        a_low  = qdd_getlow(node_a);
        a_high = qdd_gethigh(node_a);
        // Pass edge weight of current edge down to low/high
        AMP new_amp_low  = Cmul(amp_a, QDD_AMP(a_low));
        AMP new_amp_high = Cmul(amp_a, QDD_AMP(a_high));
        a_low  = qdd_bundle_ptr_amp(QDD_PTR(a_low),  new_amp_low);
        a_high = qdd_bundle_ptr_amp(QDD_PTR(a_high), new_amp_high);
    }
    if(QDD_PTR(b) == QDD_TERMINAL){
        var_b = 9000; // TODO: do this slightly cleaner
    }
    else{
        qddnode_t node_b = QDD_GETNODE(QDD_PTR(b));
        var_b  = qdd_getvar(node_b);
        b_low  = qdd_getlow(node_b);
        b_high = qdd_gethigh(node_b);
        // Pass edge weight of current edge down to low/high
        AMP new_amp_low  = Cmul(amp_b, QDD_AMP(b_low));
        AMP new_amp_high = Cmul(amp_b, QDD_AMP(b_high));
        b_low  = qdd_bundle_ptr_amp(QDD_PTR(b_low),  new_amp_low);
        b_high = qdd_bundle_ptr_amp(QDD_PTR(b_high), new_amp_high);
    }


    // Base/terminal case: same target and same variable
    if(QDD_PTR(a) == QDD_PTR(b) && var_a == var_b){
        AMP sum = Cadd(amp_a, amp_b);
        QDD res = qdd_bundle_ptr_amp(QDD_PTR(a), sum);
        return res;
    }



    sylvan_stats_count(QDD_PLUS);

    bool cachenow = 1;
    if (cachenow) {
        QDD res;
        if (cache_get3(CACHE_QDD_PLUS, sylvan_false, a, b, &res)) {
            //printf("\nlooked something up instead of recomputing for PLUS\n\n");
            sylvan_stats_count(QDD_PLUS_CACHED);
            return res;
        }
    }
    
    // Recursive case
    // If var_a != var_b, the higher one "skips" the lower variable and 
    // is effectively a don't-care for that lower variable. In this case
    // we (effectively) "reinsert" this don't-care.
    BDDVAR level = var_a;
    if(var_a > var_b){
        a_low  = a;
        a_high = a;
        level  = var_b;
    }
    else if(var_a < var_b){
        b_low  = b;
        b_high = b;
        level  = var_a;
    }

    QDD res_low, res_high;

    bdd_refs_spawn(SPAWN(qdd_plus, a_high, b_high));
    res_low = CALL(qdd_plus, a_low, b_low);
    bdd_refs_push(res_low);
    res_high = bdd_refs_sync(SYNC(qdd_plus));
    bdd_refs_pop(1);


    QDD res = qdd_makenode(level, res_low, res_high);

    if (cachenow) {
        if (cache_put3(CACHE_QDD_PLUS, sylvan_false, a, b, res)) sylvan_stats_count(QDD_PLUS_CACHEDPUT);
    }

    return res;
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

    AMP a = C_ONE;
    for (;;) {
        a = Cmul(a, QDD_AMP(q));
        
        // if the current edge is pointing to the terminal, we're done.
        if (QDD_PTR(q) == QDD_TERMINAL) break;

        // now we need to choose low or high edge of next node
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        BDDVAR var     = qdd_getvar(node);

        // Condition low/high choice on basis state vector[var]
        if (basis_state[var] == 0)
            q = node->low;
        else
            q = node->high;
    }

    // TODO: return complex struct instead of the index?
    return a;
}

QDD
create_all_zero_state(int n)
{
    assert(n >= 1);

    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    return create_basis_state(n, x);
}

QDD
create_basis_state(int n, bool* x)
{
    assert(n >= 1);

    struct qddnode node;
    PTR low_child, high_child;
    AMP low_amp, high_amp;

    // start at terminal, and build backwards
    PTR prev = QDD_TERMINAL;
    for(int k = n-1; k >=0; k--){

        if(x[k] == 0){
            low_child  = prev;
            low_amp    = C_ONE;
            high_child = QDD_TERMINAL;
            high_amp   = C_ZERO;
        }
        else if(x[k] == 1){
            low_child  = QDD_TERMINAL;
            low_amp    = C_ZERO;
            high_child = prev;
            high_amp   = C_ONE;
        }

        // pack info into node (TODO: rename this function)
        qddnode_make(&node, k, low_child, high_child, low_amp, high_amp);

        // actually make the node (i.e. add to nodetable)
        prev = QDD_PTR(qdd_makenode(k, node.low, node.high));
    }

    QDD root_edge = qdd_bundle_ptr_amp(prev, C_ONE);
    return root_edge;
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
            if(!CexactEqual(Cvalue(amp_a), Cvalue(amp_b))){
                if(verbose){
                    _print_bitstring(x, n);
                    printf(", amp a ="); Cprint(Cvalue(amp_a));
                    printf(" != amp b ="); Cprint(Cvalue(amp_b));
                    printf("\n");
                }
                return false;
            }
        }
        else{
            if(!CapproxEqual(Cvalue(amp_a), Cvalue(amp_b))){
                if(verbose){
                    _print_bitstring(x, n);
                    printf(", amp a ="); Cprint(Cvalue(amp_a));
                    printf(" !~= amp b ="); Cprint(Cvalue(amp_b));
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
_print_bitstring(bool *x, int n)
{
    for(int k=n-1; k>=0; k--) printf("%d", x[k]);
}

void
_print_qdd(QDD q)
{
    if(QDD_PTR(q) != QDD_TERMINAL){
        qddnode_t node = QDD_GETNODE(QDD_PTR(q));
        printf("%lx\t", QDD_PTR(q));
        pprint_qddnode(node);
        _print_qdd(qdd_getlow(node));
        _print_qdd(qdd_gethigh(node));
    }
}

void
print_qdd(QDD q)
{
    printf("root edge: %lx, %lx = ",QDD_PTR(q), QDD_AMP(q));
    Cprint(Cvalue(QDD_AMP(q)));
    printf("\n");
    _print_qdd(q);
}


