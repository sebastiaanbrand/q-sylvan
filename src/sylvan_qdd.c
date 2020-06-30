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
 *  [    24 bits blank   | 8 bit c | 8 bit t1 | 8 bit t2 |  16 bits gateid    ]
 *  Set control = 0 for single qubit gates
 *  (maybe put this elsewhere?)
 */
static inline uint64_t
GATE_OPID(uint32_t gateid, BDDVAR c, BDDVAR t1, BDDVAR t2)
{
    // I don't remember why we're doing this with the inputs
    uint64_t _c = c;
    uint64_t _t1 = t1;
    uint64_t _t2 = t2;
    uint64_t res = _c<<32 | _t1<<24 | _t2<<16 | gateid;
    return res;
}


/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
qddnode_getvar(qddnode_t n)
{
    return (BDDVAR) ( (n->low >> 60) | ((n->high >> 56) & 0xf0) ); // 8 bits
}

/**
 * Gets the low edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qddnode_getlow(qddnode_t n)
{
    return (QDD) n->low & 0x07ffffffffffffff; // 59 bits
}

/**
 * Gets the high edge of `n` with the AMP and PTR information, but without the
 * (halved) variable information.
 */
static inline QDD
qddnode_gethigh(qddnode_t n)
{
    return (QDD) n->high & 0x07ffffffffffffff; // 59 bits
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
 * Gets only the AMP of the low edge of `n`.
 */
static inline AMP
qddnode_getamplow(qddnode_t n)
{
    return (AMP) QDD_AMP(n->low);
}

/**
 * Gets only the AMP of the high edge of `n`.
 */
static inline AMP
qddnode_getamphigh(qddnode_t n)
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
static void qddnode_pprint(qddnode_t n)
{
    AMP amp_low  = qddnode_getamplow(n);
    AMP amp_high = qddnode_getamphigh(n);
    printf("[var=%d, low=%lx, high=%lx, ", 
             qddnode_getvar(n),
             qddnode_getptrlow(n),
             qddnode_getptrhigh(n));
    if(amp_low == C_ZERO)      printf("a=C_ZERO, ");
    else if(amp_low == C_ONE)  printf("a=C_ONE, ");
    else {
        printf("a=%lx, ",amp_low);
        printf("("); Cprint(Cvalue(amp_high)); printf(")");
    }                      
    if(amp_high == C_ZERO)     printf("b=C_ZERO ");
    else if(amp_high == C_ONE) printf("b=C_ONE, ");
    else {                     
        printf("b=%lx", amp_high);
        printf("("); Cprint(Cvalue(amp_high)); printf(")");
    }
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

static QDD // (PTR and AMP, but the amp is the norm weight from below)
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
    qdd_unmark_rec(qddnode_getlow(n));
    qdd_unmark_rec(qddnode_gethigh(n));
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
    return 1 + qdd_nodecount_mark(qddnode_getlow(n)) + qdd_nodecount_mark(qddnode_gethigh(n));
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
    BDDVAR var = qddnode_getvar(node);
    
    // "above" the desired qubit in the QDD
    if(var < qubit){

        bool cachenow = 1;
        if (cachenow) {
            QDD res;
            // check if this calculation has already been done before for this node/gate
            if (cache_get3(CACHE_QDD_GATE, GATE_OPID(gate, 0, qubit, 0), q, sylvan_false, &res)) {
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
            if (cache_put3(CACHE_QDD_GATE, GATE_OPID(gate, 0, qubit, 0), q, sylvan_false, res)) sylvan_stats_count(QDD_GATE_CACHEDPUT);
        }
        return res;
    }
    else {
        PTR l, h;
        AMP a, b;
        // exactly at qubit
        if(var == qubit){
            l = qddnode_getptrlow(node);
            h = qddnode_getptrhigh(node);
            a = qddnode_getamplow(node);
            b = qddnode_getamphigh(node);
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
    BDDVAR var = qddnode_getvar(node);

    // "above" the desired qubit in the QDD (this is where the recursive stuff
    // of cgate happens, once the control qubit has been reached, the recursive
    // search for target qubit continues in qdd_gate function)
    if(var < c){
        QDD res_low, res_high;

        bool cachenow = 1;
        if (cachenow) {
            QDD res;
            // check if this calculation has already been done before for this node/gate
            if (cache_get3(CACHE_QDD_CGATE, GATE_OPID(gate, c, t, 0), q, sylvan_false, &res)) {
                sylvan_stats_count(QDD_CGATE_CACHED);
                return res;
            }
        }

        bdd_refs_spawn(SPAWN(qdd_cgate, node->high, gate, c, t));
        res_low = CALL(qdd_cgate, node->low, gate, c, t);
        bdd_refs_push(res_low);
        res_high = bdd_refs_sync(SYNC(qdd_cgate));
        bdd_refs_pop(1);

        QDD res  = qdd_makenode(var, res_low, res_high);
        AMP new_root_amp = Cmul(QDD_AMP(q), QDD_AMP(res));
        res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);

        if (cachenow) {
            if (cache_put3(CACHE_QDD_CGATE, GATE_OPID(gate, c, t, 0), q, sylvan_false, res)) sylvan_stats_count(QDD_CGATE_CACHEDPUT);
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
        var_a  = qddnode_getvar(node_a);
        a_low  = qddnode_getlow(node_a);
        a_high = qddnode_gethigh(node_a);
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
        var_b  = qddnode_getvar(node_b);
        b_low  = qddnode_getlow(node_b);
        b_high = qddnode_gethigh(node_b);
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

QDD
qdd_swap_gate(QDD qdd, BDDVAR qubit1, BDDVAR qubit2)
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

TASK_IMPL_4(QDD, qdd_cswap_gate, QDD, qdd, BDDVAR, c, BDDVAR, t1, BDDVAR, t2)
{
    assert (c  < t1);
    assert (t1 < t2);

    bool cachenow = 1; 
    if (cachenow) {
        QDD res;
        if (cache_get3(CACHE_QDD_CGATE, GATE_OPID(GATEID_swap, c, t1, t2), qdd, 0, &res)) {
            sylvan_stats_count(QDD_CGATE_CACHED);
            return res;
        }
    }

    // similar to normal control gate, TODO: maybe generalize qdd_cgate?
    bool skipped = false;
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
        if (qddnode_getvar(node) > c) {
            skipped = true;
        }
    }

    BDDVAR var;
    QDD low, high;
    bool control_here = false;
    if (skipped) { // var > control
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        var  = c;
        control_here = true;
    }
    else {
        // not skipped, either (var == control) or (var < control)
        qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
        var  = qddnode_getvar(node);
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
        if (var == c) control_here = true;
    }

    if (control_here) {
        // apply swap to high, but not to low
        high = qdd_swap_gate(high, t1, t2);
    }
    else {
        // recursive call to both children
        bdd_refs_spawn(SPAWN(qdd_cswap_gate, high, c, t1, t2));
        low = CALL(qdd_cswap_gate, low, c, t1, t2);
        bdd_refs_push(low);
        high = bdd_refs_sync(SYNC(qdd_cswap_gate));
        bdd_refs_pop(1);
    }

    QDD res = qdd_makenode(var, low, high); 

    if (cachenow) {
        if (cache_put3(CACHE_QDD_CGATE, GATE_OPID(GATEID_swap, c, t1, t2), qdd, 0, res))
            sylvan_stats_count(QDD_CGATE_CACHEDPUT);
    }
    
    AMP new_root_amp = Cmul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}


TASK_IMPL_4(QDD, qdd_all_control_phase, QDD, qdd, BDDVAR, k, BDDVAR, n, bool*, x)
{
    // TODO: remove LACE, no branching in this function
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
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }

    // terminal case, apply phase depending on x[k] (control k on 0 or 1)
    if (k == (n-1)) {
        if (x[k] == 1) {
            AMP new_amp = Cmul(QDD_AMP(high), Clookup(Cmake(-1.0, 0.0)));
            high = qdd_bundle_ptr_amp(QDD_PTR(high), new_amp);
        }
        else {
            AMP new_amp = Cmul(QDD_AMP(low), Clookup(Cmake(-1.0, 0.0)));
            low = qdd_bundle_ptr_amp(QDD_PTR(low), new_amp);
        }
    }
    // non terminal case, choose low/high depending on x[k] (control k on 0 or 1)
    else {
        if (x[k] == 1) {
            k++; // next level
            high = CALL(qdd_all_control_phase, high, k, n, x);
            k--;
        }
        else {
            k++;
            low = CALL(qdd_all_control_phase, low, k, n, x);
            k--;
        }
    }

    QDD res = qdd_makenode(k, low, high);

    // multiply by existing edge weight on qdd
    AMP new_root_amp = Cmul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}


QDD
qdd_QFT(QDD qdd, BDDVAR first, BDDVAR last)
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

    // swap qubit order
    int num_qubits = (last - first) + 1;
    for (int j = 0; j < (int)(num_qubits/2); j++) {
        a = first + j;
        b = last  - j;
        res = qdd_swap_gate(res, a, b);
    }
    return res;
}

QDD
qdd_QFT_inv(QDD qdd, BDDVAR first, BDDVAR last)
{
    LACE_ME;
    
    int k;
    QDD res = qdd;
    BDDVAR a, b;

    // swap gates
    int num_qubits = (last - first) + 1;
    for (int j = 0; j < (int)(num_qubits/2); j++) {
        a = first + j;
        b = last  - j;
        res = qdd_swap_gate(res, a, b);
    }
    
    // H gates and phase gates (but now backwards)
    for (a = last + 1; a-- > 0; ) { // weird for loop because BDDVARs are unsigned

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
_qdd_grover_iteration(QDD qdd, BDDVAR n, bool* flag)
{
    LACE_ME; // TODO: change stuff so we don't need to LACE_ME in every function

    // "oracle" call  (apply -1 flag to desired amplitude)
    qdd = qdd_all_control_phase(qdd, n, flag);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Phase on all amplitudes except |000...0>
    bool x[n]; 
    for(BDDVAR k = 0; k < n; k++) x[k] = 0;
    qdd = qdd_all_control_phase(qdd, n, x);
    AMP new_root_amp = Cmul(QDD_AMP(qdd), Clookup(Cmake(-1.0, 0.0)));
    qdd = qdd_bundle_ptr_amp(QDD_PTR(qdd), new_root_amp);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    return qdd;
}

QDD
qdd_grover(BDDVAR n, bool* flag)
{   
    LACE_ME;

    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // start with all zero state |000...0>
    QDD qdd = qdd_create_all_zero_state(n);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    for (uint32_t i = 0; i < R; i++) qdd = _qdd_grover_iteration(qdd, n, flag);

    return qdd;
}

/*******************************<Shor components>******************************/

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
        for (int j = i; j <= num_qubits; j++){
            if (a[j] == 1) {
                k = (j - i) + 1;
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
        for (int j = i; j <= num_qubits; j++){
            if (a[j] == 1) {
                k = (j - i) + 1;
                res = qdd_gate(res, GATEID_Rk_dag(k), qubit);
            }
        }
    }
    return res;
}



/******************************</Shor components>******************************/

QDD
qdd_measure_q0(QDD qdd, int *m, double *p)
{  
    LACE_ME;

    // get probabilities for q0 = |0> and q0 = |1>
    qddnode_t node;
    bool skipped = false;
    if (QDD_PTR(qdd) == QDD_TERMINAL) {
        skipped = true;
    }
    else {
        node = QDD_GETNODE(QDD_PTR(qdd));
        if (qddnode_getvar(node) != 0) {
            skipped = true;
        }
    }
    QDD low, high;
    if (skipped) {
        // if skipped q0 is a don't care, treat separately?
        low  = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
        high = qdd_bundle_ptr_amp(QDD_PTR(qdd), C_ONE);
    }
    else {
        low  = qddnode_getlow(node);
        high = qddnode_gethigh(node);
    }
    double prob_low  = _qdd_unnormed_prob(low);
    double prob_high = _qdd_unnormed_prob(high);
    double prob_root = _prob(QDD_AMP(qdd));
    prob_low  *= prob_root; // printf("p low  = %.60f\n", prob_low);
    prob_high *= prob_root; // printf("p high = %.60f\n", prob_high);
    assert(abs(prob_low + prob_high - 1.0) < TOLERANCE);

    // flip a coin
    float rnd = ((float)rand())/RAND_MAX;
    *m = (rnd < prob_low) ? 0 : 1;
    *p = prob_low;

    // produce post-measurement state
    AMP norm, normalized;
    if (*m == 0) {
        high       = QDD_TERMINAL;
        norm       = Clookup(Cmake(sqrt(prob_low), 0.0));
        normalized = Cdiv(QDD_AMP(low), norm);
        low        = qdd_bundle_ptr_amp(QDD_PTR(low), normalized);
    }
    else {
        low        = QDD_TERMINAL;
        norm       = Clookup(Cmake(sqrt(prob_high), 0.0));
        normalized = Cdiv(QDD_AMP(high), norm);
        high       = qdd_bundle_ptr_amp(QDD_PTR(high), normalized);
    }

    QDD res = qdd_makenode(0, low, high);
    AMP new_root_amp = Cmul(QDD_AMP(qdd), QDD_AMP(res));
    res = qdd_bundle_ptr_amp(QDD_PTR(res), new_root_amp);
    return res;
}


TASK_IMPL_1(QDD, _qdd_unnormed_prob, QDD, qdd)
{
    if (QDD_PTR(qdd) == QDD_TERMINAL) return _prob(QDD_AMP(qdd));

    bool cachenow = 1;
    if (cachenow) {
        QDD res;
        // check if this calculation has already been done before for this node/gate
        if (cache_get3(CACHE_QDD_PROB, 0LL, qdd, 0LL, &res)) {
            sylvan_stats_count(QDD_GATE_CACHED);
            return res;
        }
    }
    
    qddnode_t node = QDD_GETNODE(QDD_PTR(qdd));
    double p_low, p_high, res;

    bdd_refs_spawn(SPAWN(_qdd_unnormed_prob, qddnode_gethigh(node)));
    p_low = CALL(_qdd_unnormed_prob, qddnode_getlow(node));
    bdd_refs_push(p_low); // Q: this is not a bdd/qdd node, do we need to protect this?
    p_high = bdd_refs_sync(SYNC(_qdd_unnormed_prob)); // syncs SPAWN
    bdd_refs_pop(1);

    res = (p_low + p_high) * _prob(QDD_AMP(qdd));

    if (cachenow) {
        if (cache_put3(CACHE_QDD_PROB, 0LL, qdd, 0LL, res)) sylvan_stats_count(QDD_PROB_CACHEDPUT);
    }

    return res;
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
        BDDVAR var     = qddnode_getvar(node);

        // Condition low/high choice on basis state vector[var]
        if (basis_state[var] == 0)
            q = node->low;
        else
            q = node->high;
    }

    // TODO: return complex struct instead of the index?
    return a;
}

double
_prob(AMP a) 
{
    // move to qdd_int file?
    complex_t c = Cvalue(a);
    double abs = sqrt( (c.r*c.r) + (c.i*c.i) );
    return (abs*abs);
}

QDD
qdd_create_all_zero_state(int n)
{
    assert(n >= 1);

    bool x[n];
    for(int k=0; k<n; k++) x[k] = 0;
    return qdd_create_basis_state(n, x);
}

QDD
qdd_create_basis_state(int n, bool* x)
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
        if(!qddnode_getmark(node)){
            qddnode_setmark(node, 1);
            printf("%lx\t", QDD_PTR(q));
            qddnode_pprint(node);
            _print_qdd(qddnode_getlow(node));
            _print_qdd(qddnode_gethigh(node));
        }
    }
}

void
qdd_printnodes(QDD q)
{
    printf("root edge: %lx, %lx = ",QDD_PTR(q), QDD_AMP(q));
    Cprint(Cvalue(QDD_AMP(q)));
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
        complex_t val = Cvalue(a);
        if (val.r != 0.0) fprintf(out, "%.3lf", val.r);
        if (val.i > 0.0) fprintf(out, "+%.3lfi", val.i);
        else if (val.i < 0.0) fprintf(out, "-%.3lfi", val.i);
    }
    fprintf(out, "\"");
}

static void
qdd_fprintdot_rec(FILE *out, QDD qdd)
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
    qdd_fprintdot_rec(out, qddnode_getlow(n));
    qdd_fprintdot_rec(out, qddnode_gethigh(n));

    // add edge from this node to each child
    fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed",
                QDD_PTR(qdd), qddnode_getptrlow(n));
    qdd_fprintdot_label(out, qddnode_getamplow(n));
    fprintf(out, "];\n");
    fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid",
                QDD_PTR(qdd), qddnode_getptrhigh(n));
    qdd_fprintdot_label(out, qddnode_getamphigh(n));
    fprintf(out, "];\n");
    
    // TODO: edge weights
}

void
qdd_fprintdot(FILE *out, QDD qdd)
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
    qdd_fprintdot_rec(out, qdd);
    qdd_unmark_rec(qdd);

    fprintf(out, "}\n");
}


