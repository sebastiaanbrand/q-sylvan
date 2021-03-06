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

/* Do not include this file directly. Instead, include sylvan.h */


/**
 *
A QMDD can represent exponentially-sized vectors and matrices of complex numbers
including quantum states and gates. For those familiar with the literature on
decision diagrams, it might be helpful to understand the relation between QMDDDs
and other DD variants. Perhaps surprisingly, the QMDD data structure is more
similar to binary DDs (but additionally adorned with edge weights) than to a
multi-terminal or multi-values DDs, as its name and purpose seems to suggest.
A multi-terminal DD, aka algebraic DD or ADD, represents a function
f : B^n --> D, where D,  the domain, is an algebra, e.g., the complex numbers.
ADDs can therefore immediately be used to encode quantum systems as well,
though less efficiently [1].
A multi-valued DD represents functions f : D^n --> B and is typically used to
encode high-level systems such as programs with multi-valued variables.

Structurally, a QMDD is an Ordered and Reduced BDD (ORBDD) with weighted edges.
The root node also has a weighted start edge. As a BDD, each node has a
variable x and two outgoing edges to QMDDs: a low and a high edge representing
x=0 and x=1. Edge weights are normalized by dividing the weights of the low and
high edges by that of the low edge (unless it equals zero) and propagating alpha
upwards through the incoming edge. This normalization achieves additional
reductions in the structure and should not be confused with the normalization
of quantum states, i.e., QMDDs are more general than quantum systems.

Formally, a QMDD edge e is a pair (Q,w), where Q is (a pointer to) a QMMD node,
w in C is a weight representing a normalized amplitude. A QMDD node is either a
constant 1 or a tuple (x,l,h) where x is a (Boolean) variable and
l, h QMDD edges such that:
- x < l.Q.x, r.Q.x,         (similar to the ordering in ORBDDs)
- (x,l,h) is unique,        (similar to the reduction in ORBDDs)
- l != h, and               (similar to the deletion rule from ORBDDs)
- l.w in {0,1}.             (by normalization)

A QMDD (edge) e is thus interpreted as a function [[ e ]] : (X --> B) --> C,
where X are the variables, and B,C the Boolean / Complex domains.
Each path to the 1 leaf node represents an assignment A : X --> B to the
variables x in X of the nodes along the path (x=0 for low and x=1 for high),
and an amplitude defined by multiplying the weights along the path.
Formally, the semantics of [[ e ]] for a QMDD edge e is defined as:
[[ (1,w) ]]      := { {} --> w }
[[ (Q,w) ]]      := { {Q.x = 0} U A --> w * Q.l.w | (s,A) in [[ Q.l ]]  } U
                    { {Q.x = 1} U A --> w * Q.r.w | (s,A) in [[ Q.h ]]  }
For each assignment A in B^X, we may interpret [[ e ]](A) as usual by
interpreting unassigned variables as `don't cares`, i.e.,
 [[ e ]](A)      := w such that exists (A',w) in [[ e ]] with A' subseteq A.
 (Note that by the ordered property, such an (A',w) is unqiue in [[ e ]].)
Alternative, we might extend the QMDD of e by re-adding all deleted nodes (see
definition above). This can be done by adding `don't care` nodes (x,(Q,1),(Q,1))
between edges which skip a variable, i.e., l.x != x+1 or h.x != x+1.
Formally, we replace each node (x,l,h) with (x, n_x+1, h), where
- n_i := (i, (1,n_i+1),  (1,n_i+1) )            for i < l.x - 1, and
- n_i := (i, l, h)                              for i = l.x - 1.
And analogously for the high edges.

 */


#ifndef SYLVAN_QDD_H
#define SYLVAN_QDD_H

#include <stdbool.h>
#include <stdint.h>

#include "amp_storage/flt.h"
#include "sylvan_qdd_gates.h"
#include "sylvan_mtbdd.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// TODO: organize this file better
/**
 * 
 * QDD : essentially an edge, containing [var, amp, pointer]
 *  - var: the variable/qubit numer (unclear whether this variable refers to the
 *         node this edge is an outgoing edge of, or the node this edge is an 
 *         incomming edge to.)
 *  - amp: the (index of) the amplitude on this edge.
 *  - pointer: a number which identifies the node this edge is pointing _to_.
 * qddnote_t : essentially a node, containing two edges (to its children).
 */

typedef uint64_t QDD; // QDD edge (contains AMP and PTR)
typedef uint64_t AMP; // amplitude index
typedef uint64_t PTR; // node index

typedef struct __attribute__((packed)) qddnode {
    QDD low, high;
} *qddnode_t; // 16 bytes

static const PTR        QDD_TERMINAL = 1;
static const BDDVAR     QDD_INVALID_VAR = UINT8_MAX;

typedef enum weight_norm_strategy {
    NORM_LOW,
    NORM_LARGEST,
    NORM_SUM,
    n_norm_stragegies
} weight_norm_strategy_t;

/**
 * Similar initialization as for MTBDDs + amplitude table init.
 * Setting tolerance to -1 uses default tolerance.
 * real table: stores 2 real values per edge weight, instead of 1 tuple
 * NOTE: this function doesn't currently check if the combination of table
 * sizes (edge weight table + node table) works in combination with using 
 * a real-table or complex-table.
 */
void sylvan_init_qdd(size_t ctable_size, double ctable_tolerance, int amps_backend, int norm_strat);
void sylvan_init_qdd_defaults(size_t ctable_size);
void qdd_set_testing_mode(bool on);
void qdd_set_caching_granularity(int granularity);


/*******************<garbage collection, references, marking>******************/

VOID_TASK_DECL_1(qdd_gc_mark_rec, QDD);
#define qdd_gc_mark_rec(qdd) CALL(qdd_gc_mark_rec, qdd)

/**
 * Store the pointer <ptr> in the pointers table.
 */
void qdd_protect(QDD* ptr);

/**
 * Delete the pointer <ptr> from the pointers table.
 */
void qdd_unprotect(QDD* ptr);

/**
 * Compute the number of pointers in the pointers table.
 */
size_t qdd_count_protected(void);

/**
 * Push a QDD variable to the pointer reference stack.
 * During gc the variable will be inspected and the contents will be marked.
 */
void qdd_refs_pushptr(const QDD *ptr);

/**
 * Pop the last <amount> QDD variables from the pointer reference stack.
 */
void qdd_refs_popptr(size_t amount);

/**
 * Push a QDD to the values reference stack.
 * During garbage collection the references QDD will be marked.
 */
QDD qdd_refs_push(QDD qdd);

/**
 * Pop the last <amount> QDDs from the values reference stack.
 */
void qdd_refs_pop(long amount);

/**
 * Push a Task that returns a QDD to the tasks reference stack.
 * Usage: qdd_refs_spawn(SPAWN(function, ...));
 */
void qdd_refs_spawn(Task *t);

/**
 * Pop a Task from the task reference stack.
 * Usage: QDD result = qdd_refs_sync(SYNC(function));
 */
QDD qdd_refs_sync(QDD qdd);

/******************</garbage collection, references, marking>******************/


/*******************************<applying gates>*******************************/

/* Applies given (single qubit) gate to |q>. (Wrapper function) */
#define qdd_gate(qdd,gate,target) (CALL(qdd_gate,qdd,gate,target))
TASK_DECL_3(QDD, qdd_gate, QDD, gate_id_t, BDDVAR);

/* Applies given controlled gate to |q>. (Wrapper function) */
#define qdd_cgate(qdd,gate,c,t) (CALL(qdd_cgate,qdd,gate,c,t))
#define qdd_cgate2(qdd,gate,c1,c2,t) (CALL(qdd_cgate2,qdd,gate,c1,c2,t))
#define qdd_cgate3(qdd,gate,c1,c2,c3,t) (CALL(qdd_cgate3,qdd,gate,c1,c2,c3,t))
TASK_DECL_4(QDD, qdd_cgate,  QDD, gate_id_t, BDDVAR, BDDVAR);
TASK_DECL_5(QDD, qdd_cgate2, QDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR);
TASK_DECL_6(QDD, qdd_cgate3, QDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR, BDDVAR);

/* Applies given controlled gate to |q>. (Wrapper function) */
#define qdd_cgate_range(qdd,gate,c_first,c_last,t) (CALL(qdd_cgate_range,qdd,gate,c_first,c_last,t))
TASK_DECL_5(QDD, qdd_cgate_range, QDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR);


/**
 * Recursive implementation of vector addition.
 */
#define qdd_plus(a,b) (CALL(qdd_plus,a,b))
TASK_DECL_2(QDD, qdd_plus, QDD, QDD);

/**
 * Recursive implementation of applying single qubit gates
 */
#define qdd_gate_rec(q,gate,target) (CALL(qdd_gate_rec,q,gate,target))
TASK_DECL_3(QDD, qdd_gate_rec, QDD, gate_id_t, BDDVAR);

/**
 * Recursive implementation of applying controlled gates
 */
#define qdd_cgate_rec(q,gate,cs,t) (CALL(qdd_cgate_rec,q,gate,cs,0,t))
TASK_DECL_5(QDD, qdd_cgate_rec, QDD, gate_id_t, BDDVAR*, uint32_t, BDDVAR);

/**
 * Recursive implementation of applying controlled gates where the controlles 
 * are defined by a range 'c_first' through 'c_last'.
 */
#define qdd_cgate_range_rec(q,gate,c_first,c_last,t) (CALL(qdd_cgate_range_rec,q,gate,c_first,c_last,t,0))
TASK_DECL_6(QDD, qdd_cgate_range_rec, QDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR, BDDVAR);


/* Computes Mat * |vec> (Wrapper function) */
#define qdd_matvec_mult(mat,vec,nvars) (CALL(qdd_matvec_mult,mat,vec,nvars))
TASK_DECL_3(QDD, qdd_matvec_mult, QDD, QDD, BDDVAR);

/* Computes A*B (note generally AB != BA) (Wrapper function) */
#define qdd_matmat_mult(a,b,nvars) (CALL(qdd_matmat_mult,a,b,nvars))
TASK_DECL_3(QDD, qdd_matmat_mult, QDD, QDD, BDDVAR);

/**
 * Recursive implementation of matrix-vector mult and matrix-matrix mult.
 */
TASK_DECL_4(QDD, qdd_matvec_mult_rec, QDD, QDD, BDDVAR, BDDVAR);
TASK_DECL_4(QDD, qdd_matmat_mult_rec, QDD, QDD, BDDVAR, BDDVAR);

/* Multiply some qdd by a scalar. */
QDD qdd_scalar_mult(QDD qdd, complex_t c);

/**
 * Increases all the variable number in QDD a by k (used for tensor product)
 * 
 * @param a QDD over n vars {j, j+1, ..., j+n-1} (generally j=0)
 * 
 * @return QDD over n vars {j+k, j+k+1, ..., j+k+n-1}
 */
QDD qdd_increase_all_vars(QDD a, int k);

/* Replace the terminal node in a with b (effectively stacks a and b) 
* (used for tensor product)
*/
QDD qdd_replace_terminal(QDD a, PTR b);

/**
 * @param a QDD over vars 0...n-1 (n = nvars_a)
 * @param b QDD over vars 0...m-1
 * @param nvars_a number of vars of QDD a
 * 
 * @return QDD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
QDD qdd_tensor_prod(QDD a, QDD b, BDDVAR nvars_a);

/**
 * Computes the tensor product of vec (tensor) vec
 * 
 * @param a QDD over vars 0...n-1 (n = nvars_a)
 * @param b QDD over vars 0...m-1
 * @param nvars_a number of vars of QDD a
 * 
 * @return QDD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
#define qdd_vec_tensor_prod(a, b, nvars_a) qdd_tensor_prod(a,b,nvars_a)

/**
 * Computes the tensor product of mat (tensor) mat
 * 
 * @param a QDD over vars 0...2n-1 (n = nvars_a)
 * @param b QDD over vars 0...2m-1
 * @param nvars_a number of vars of QDD a (counting only unprimed)
 * 
 * @return QDD over vars 0...2n-1...(2n+2m)-1, representing a (tensor) b
 */
#define qdd_mat_tensor_prod(a, b, nvars_a) qdd_tensor_prod(a,b,2*nvars_a)

/******************************</applying gates>*******************************/



/*********************<applying (controlled) sub-circuits>*********************/

// Circuit IDs
#define CIRCID_swap          0
#define CIRCID_reverse_range 1
#define CIRCID_QFT           2
#define CIRCID_QFT_inv       3

// For now we have at most 3 control qubits
static const uint32_t MAX_CONTROLS = 3;

/**
 * Circuit which implements a SWAP gate from single-qubit and controlled gates.
 */
QDD qdd_circuit_swap(QDD qdd, BDDVAR qubit1, BDDVAR qubit2);

/**
 * Circuit which reverses the order of the qubits in the given range.
 */
QDD qdd_circuit_reverse_range(QDD qdd, BDDVAR first, BDDVAR last);

/**
 * Executes the QFT circuit on qubits `first` through `last`.
 */
QDD qdd_circuit_QFT(QDD qdd, BDDVAR first, BDDVAR last);

/**
 * Executes the inverse QFT circuit on qubits `first` through `last`.
 */
QDD qdd_circuit_QFT_inv(QDD qdd, BDDVAR first, BDDVAR last);

/**
 * Applies the given circuit (parameters can be two targets or a range
 * depending on the circuit.)
 */
QDD qdd_circuit(QDD qdd, uint32_t circ_id, BDDVAR t1, BDDVAR t2);

/**
 * Generalized implementation of applying controlled versions of sub-circuit
 * functions defined here.
 * @param circ_id CIRCID_something
 * @param cs BDDVAR[] of control qubits. Needs to be length 3. If using fewer
 *           controls use e.g. cs = [c1, c2, QDD_INVALID_VAR]
 * @param t1 BDDVAR. Parameter 1 for given circuit.
 * @param t2 BDDVAR. Parameter 2 for given circuit.
 */
#define qdd_ccircuit(qdd, circ_id, cs, t1, t2) (CALL(qdd_ccircuit,qdd,circ_id,cs,0,t1,t2));
TASK_DECL_6(QDD, qdd_ccircuit, QDD, uint32_t, BDDVAR*, uint32_t, BDDVAR, BDDVAR);

/**
 * Applies a phase of -1 to a single basis state |x>.
 * This is a CZ gate where we control on all qubits and when x_k = 0 we control
 * qubit k on 0, and where x_k = 1 we control qubit k on 1.
 * 
 * @param qdd A QDD encoding some quantum state |\psi>.
 * @param n Number of qubits.
 * @param x A bitstring x of some computational basis state |x>.
 * 
 * TODO: generalize this to control on some but not all qubits.
 */
QDD qdd_all_control_phase(QDD qdd, BDDVAR n, bool *x);

/********************</applying (controlled) sub-circuits>*********************/


/**
 * Computational basis measurement on qubit q_k.
 * 
 * @param qdd A QDD encoding of some n qubit state.
 * @param k Which qubit to measure.
 * @param m Return of measurement outcome (0 or 1).
 * @param p Return of measurement probability.
 * 
 * @return QDD of post-measurement state corresponding to measurement outcome.
 */
QDD qdd_measure_qubit(QDD qqd, BDDVAR k, BDDVAR nvars, int *m, double *p);
QDD qdd_measure_q0(QDD qdd, BDDVAR nvars, int *m, double *p);

/**
 * Computational basis measurement of all n qubits in the qdd.
 * 
 * @param qdd A QDD encoding an n qubit state |\psi>.
 * @param n Number of qubits.
 * @param ms Array of lenght n where the measurement outcomes are put.
 * @param p Return of measurement probability |<\psi|ms>|^2.
 * 
 * @return QDD of post-measurement state (computational basis state |ms>).
 */
QDD qdd_measure_all(QDD qdd, BDDVAR n, bool* ms, double *p);

/**
 * (Recursive) helper function for obtaining probabilities for measurements
 */
#define qdd_unnormed_prob(qdd, topvar, nvars) (CALL(qdd_unnormed_prob,qdd,topvar,nvars));
TASK_DECL_3(double, qdd_unnormed_prob, QDD, BDDVAR, BDDVAR);

/**
 * Get the index of the amplitude of given basis state.
 * 
 * @param qdd A QDD encoding some quantum state |\psi>.
 * @param basis_state A bitstring x of some computational basis state |x>.
 * 
 * @return The index of amplitude <x|\psi>.
 */
AMP qdd_get_amplitude(QDD qdd, bool* basis_state);

/**
 * Get amplitude of given basis state.
 * 
 * @param qdd A QDD encoding some quantum state |\psi>.
 * @param basis_state A bitstring x of some computational basis state |x>.
 * 
 * @return The amplitude <x|\psi>.
 */
complex_t qdd_get_amplitude_as_complex(QDD qdd, bool *basis_state);

/**
 * Creates a QDD for an n-qubit state |00...0>.
 * 
 * @param n Number of qubits.
 * 
 * @return A QDD encoding the n-qubit state |00..0>.
 */
QDD qdd_create_all_zero_state(BDDVAR n);

/**
 * Creates a QDD for an n-qubit state |x>.
 * 
 * @param n Number of qubits.
 * @param x A bitstring x \in {0,1}^n. 
 * 
 * @return A QDD encoding of the n-qubit state |x>.
 */
QDD qdd_create_basis_state(BDDVAR n, bool* x);

/**
 * Creates a QDD representing the matrix I \tensor I \tensor ... \tensor I.
 * 
 * @param n Number of qubits.
 * 
 * @return A QDD encoding of I \tensor I \tensor ... \tensor I
 */
QDD qdd_create_all_identity_matrix(BDDVAR n);


/**
 * Creates a QDD matrix which applies gate U to qubit t and I to all others.
 * 
 * @param n Total number of qubits.
 * @param t Target qubit for given gate.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A QDD encoding of I_0 \tensor ... U_t ... \tensor I_{n-1}
 */
QDD qdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid);

/**
 * Creates a QDD matrix which applies the given list of n gates to n qubits.
 * 
 * @param n Total number of qubits.
 * @param gateids List of length n of gate ID of predefined single qubit gates.
 * 
 * @return A QDD encoding of U_0 \tensor U_1 \tensor ... \tensor U_{n-1}
 */
QDD qdd_create_single_qubit_gates(BDDVAR n, gate_id_t *gateid);

/**
 * Creates a QDD matrix which applies single qubit gate U to all qubits.
 * 
 * @param n Total number of qubits.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A QDD encoding of U^{\tensor n}
 */
QDD qdd_create_single_qubit_gates_same(BDDVAR n, gate_id_t gateid);

/**
 * Creates a QDD matrix which does gateid(c,t) and I on all other qubits.
 * 
 * @param n Total number of qubits.
 * @param c Control qubit (for now we assume/need that c < t).
 * @param t Target qubit for given gate.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A matrix QDD encoding of the controlled gate on given qubits.
 */
QDD qdd_create_controlled_gate(BDDVAR n, BDDVAR c, BDDVAR t, gate_id_t gateid);

/**
 * Creates a controlled-`gateid` gate which acts on the qubits as specified by
 * `c_options`.
 * 
 * @param n Total number of qubits.
 * @param c_options Array of length n with option for each qubit k: {
 *        -1 : ignore qubit k (apply I), 
 *         0 : control on q_k = |0>
 *         1 : control on q_k = |1>
 *         2 : target qubit }
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A matrix QDD encoding of the multi-controlled gate on given qubits.
 */
#define qdd_create_multi_cgate(n,c_options,gateid) qdd_create_multi_cgate_rec(n,c_options,gateid,0)
QDD qdd_create_multi_cgate_rec(BDDVAR n, int *c_options, gate_id_t gateid, BDDVAR k);

/**
 * Creates an n-qubit controlled Z gate, controlled on all qubits. 
 * 
 * @param n Number of qubits.
 * @param x A bitstring x whether to control qubit k on x[k] = 0 or x[k] = 1.
 * 
 * @return A matrix QDD encoding of an n-qubit controlled Z gate.
 */
QDD qdd_create_all_control_phase(BDDVAR n, bool *x);


/**
 * Removes any global phase from state vector, assumes `qdd` is the root edge.
 */
QDD qdd_remove_global_phase(QDD qdd);

/**
 * Checks if the states represented by two QDDs are the same. In principle the
 * root edges of these QDDs should always be the same if they represent the
 * same state, so this is mostly a debug/testing function.
 * 
 * @param a A QDD representing a quantum state |psi>
 * @param b A QDD representing a quantum state |phi>
 * @param n Number of qubits
 * @param exact If true, |psi> should equal |phi> exactly (float equal),
 * otherwise allow for preset float equivalence margine.
 * 
 * @returns True iff |psi> == |phi>.
 */
bool qdd_equivalent(QDD a, QDD b, int n, bool exact, bool verbose);

/**
 * Brute force sanity check to verify a given QDD encodes (something close to) 
 * a unit vector.
 */
bool qdd_is_close_to_unitvector(QDD qdd, BDDVAR n, double tol);
bool qdd_is_unitvector(QDD qdd, BDDVAR n);

/** Sanity check to see if the QDD variables are ordered and < nvars. */
#define qdd_is_ordered(qdd,nvars) (CALL(qdd_is_ordered,qdd,0,nvars))
TASK_DECL_3(bool, qdd_is_ordered, QDD, BDDVAR, BDDVAR);

// counts the nodes by recursively marking them (and unmarks when done)
uint64_t qdd_countnodes(QDD qdd);

// temp trigger for gc of node table every n gates
void qdd_set_periodic_gc_nodetable(int every_n_gates);
/* enabled by default */
void qdd_set_auto_gc_amp_table(bool enabled);
/* default 0.5 */
void qdd_set_gc_amp_table_thres(double fraction_filled);
double qdd_get_gc_amp_table_thres();
void qdd_gc_amp_table();
void qdd_test_gc_amptable(QDD *keep);
/**
 * Recursive function for moving amps from old to new amp table.
 */
#define _fill_new_amp_table(qdd) (CALL(_fill_new_amp_table, qdd))
TASK_DECL_1(QDD, _fill_new_amp_table, QDD);

/**
 * Write a .dot representation of a given QDD
 */
void qdd_fprintdot(FILE *out, QDD qdd, bool draw_zeros);

/*******************************<logging stats>********************************/
void qdd_stats_start(FILE *out);
void qdd_stats_set_granularity(uint32_t g); // log every 'g' gates (default 1)
void qdd_stats_log(QDD qdd);
uint64_t qdd_stats_get_nodes_peak();
double qdd_stats_get_nodes_avg();
uint64_t qdd_stats_get_logcounter();
void qdd_stats_finish();
/******************************</logging stats>********************************/

// debug stuff
void qdd_printnodes(QDD q);
bool _next_bitstring(bool *x, int n);
void _print_bitstring(bool *x, int n, bool backwards);
uint64_t bitarray_to_int(bool *x, int n, bool MSB_first);
bool * int_to_bitarray(uint64_t n, int length, bool MSB_first);
bool bit_from_int(uint64_t a, uint8_t index);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
