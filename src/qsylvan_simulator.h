
/**
*
*  This code is created by Sebastiaan Brand based on the earlier developed MTBDD method by Tom van Dijk.
*
*  It is currently not used in the QSylvan top-layer, and its unit tests.
*
*/

#ifndef QSYLVAN_SIMULATOR_H
#define QSYLVAN_SIMULATOR_H

#include <sylvan_int.h>
#include <qsylvan_gates.h>

typedef AADD QMDD; // QMDD edge (contains AMP and PTR)
typedef AADD_WGT AMP; // edge weight index
typedef AADD_TARG PTR; // node index

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************************<Initialization>********************************/

/**
 * TODO: description of params.
 * 
 * This could be a constructor of the QSylvan class "QSimulator()".
 * 
*/

/*
class QSimulator()
{
    size_t weight_table_size = DEFAULT_WEIGHT_TABLE_SIZE; 
    double weight_table_tolerance = DEFAULT_TABLE_TOLERANCE; 
    int edge_weight_backend = DEFAULT_WEIGHT_BACKEND;
    int norm_strategy = DEFAULT_NORM_STRATEGY;

    QSimulator(size_t weight_table_size, double weight_table_tolerance, int edge_weight_backend, int norm_strategy)
    {
        qsylvan_init_simulator(weight_table_size, weight_table_tolerance, edge_weight_backend, norm_strategy);
        qsylvan_init_defaults(size_t wgt_tab_size);
    }

}
*/

void qsylvan_init_simulator(size_t wgt_tab_size, double wgt_tab_tolerance, int edge_weigth_backend, int norm_strat);
void qsylvan_init_defaults(size_t wgt_tab_size);

/*****************************</Initialization>********************************/


/***************************<Initial state creation>***************************/

/**
 * Creates a QMDD for an n-qubit state |00...0>.
 * 
 * @param n Number of qubits.
 * 
 * @return A QMDD encoding the n-qubit state |00..0>.
 */
QMDD qmdd_create_all_zero_state(BDDVAR n);

/**
 * Creates a QMDD for an n-qubit state |x>.
 * 
 * @param n Number of qubits.
 * @param x A bitstring x \in {0,1}^n. 
 * 
 * @return A QMDD encoding of the n-qubit state |x>.
 */
QMDD qmdd_create_basis_state(BDDVAR n, bool* x);

/**
 * Creates a QMDD representing the matrix I \tensor I \tensor ... \tensor I.
 * 
 * @param n Number of qubits.
 * 
 * @return A QMDD encoding of I \tensor I \tensor ... \tensor I
 */
QMDD qmdd_create_all_identity_matrix(BDDVAR n);

/**
 * Creates a QMDD matrix which applies gate U to qubit t and I to all others.
 * 
 * @param n Total number of qubits.
 * @param t Target qubit for given gate.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A QMDD encoding of I_0 \tensor ... U_t ... \tensor I_{n-1}
 */
QMDD qmdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid);

/**
 * Creates a QMDD matrix which applies the given list of n gates to n qubits.
 * 
 * @param n Total number of qubits.
 * @param gateids List of length n of gate ID of predefined single qubit gates.
 * 
 * @return A QMDD encoding of U_0 \tensor U_1 \tensor ... \tensor U_{n-1}
 */
QMDD qmdd_create_single_qubit_gates(BDDVAR n, gate_id_t *gateid);

/**
 * Creates a QMDD matrix which applies single qubit gate U to all qubits.
 * 
 * @param n Total number of qubits.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A QMDD encoding of U^{\tensor n}
 */
QMDD qmdd_create_single_qubit_gates_same(BDDVAR n, gate_id_t gateid);

/**
 * Creates a QMDD matrix which does gateid(c,t) and I on all other qubits.
 * 
 * @param n Total number of qubits.
 * @param c Control qubit (for now we assume/need that c < t).
 * @param t Target qubit for given gate.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return A matrix QMDD encoding of the controlled gate on given qubits.
 */
QMDD qmdd_create_controlled_gate(BDDVAR n, BDDVAR c, BDDVAR t, gate_id_t gateid);

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
 * @return A matrix QMDD encoding of the multi-controlled gate on given qubits.
 */
#define qmdd_create_multi_cgate(n,c_options,gateid) qmdd_create_multi_cgate_rec(n,c_options,gateid,0)
QMDD qmdd_create_multi_cgate_rec(BDDVAR n, int *c_options, gate_id_t gateid, BDDVAR k);

/**
 * Creates an n-qubit controlled Z gate, controlled on all qubits. 
 * 
 * @param n Number of qubits.
 * @param x A bitstring x whether to control qubit k on x[k] = 0 or x[k] = 1.
 * 
 * @return A matrix QMDD encoding of an n-qubit controlled Z gate.
 */
QMDD qmdd_create_all_control_phase(BDDVAR n, bool *x);

/**************************</Initial state creation>***************************/





/***********************<Measurements and probabilities>***********************/

/**
 * Computational basis measurement on qubit q_k.
 * 
 * @param qmdd A QMDD encoding of some n qubit state.
 * @param k Which qubit to measure.
 * @param m Return of measurement outcome (0 or 1).
 * @param p Return of measurement probability.
 * 
 * @return QMDD of post-measurement state corresponding to measurement outcome.
 */
QMDD qmdd_measure_qubit(QMDD qqd, BDDVAR k, BDDVAR nvars, int *m, double *p);
QMDD qmdd_measure_q0(QMDD qmdd, BDDVAR nvars, int *m, double *p);

/**
 * Computational basis measurement of all n qubits in the qmdd.
 * 
 * @param qmdd A QMDD encoding an n qubit state |psi>.
 * @param n Number of qubits.
 * @param ms Array of lenght n where the measurement outcomes are put.
 * @param p Return of measurement probability |<psi|ms>|^2.
 * 
 * @return QMDD of post-measurement state (computational basis state |ms>).
 */
QMDD qmdd_measure_all(QMDD qmdd, BDDVAR n, bool* ms, double *p);

/**
 * (Recursive) helper function for obtaining probabilities for measurements
 */
#define qmdd_unnormed_prob(qmdd, topvar, nvars) (RUN(qmdd_unnormed_prob,qmdd,topvar,nvars))
TASK_DECL_3(double, qmdd_unnormed_prob, QMDD, BDDVAR, BDDVAR);

/**
 * Given a state vector as QMDD, get its (L2) norm.
 */
#define qmdd_get_norm(qmdd, nvars) (RUN(qmdd_unnormed_prob,qmdd,0,nvars))

/**
 * Get amplitude of given basis state.
 * 
 * @param qmdd A QMDD encoding some quantum state |psi>.
 * @param basis_state A bitstring x of some computational basis state |x>.
 * 
 * @return The amplitude <x|psi>.
 */
complex_t qmdd_get_amplitude(QMDD qmdd, bool *basis_state);

/**
 * Computes the probability from a given edge weight index.
 * 
 * @param a Edge weight index
 * @return The probability |value(a)|^2
 */
double qmdd_amp_to_prob(AMP a);

/**
 * Turns a given probability into a quantum amplitude, and stores it in the
 * edge weight table.
 * 
 * @param a Probability
 * @return Edge weight table index of sqrt(a)
 */
AMP qmdd_amp_from_prob(double a);

/**********************</Measurements and probabilities>***********************/





/*******************************<Applying gates>*******************************/

// For now we have at most 3 control qubits
#define MAX_CONTROLS 3

/* Applies given (single qubit) gate to |q>. (Wrapper function) */
#define qmdd_gate(qmdd,gate,target) (RUN(qmdd_gate,qmdd,gate,target))
TASK_DECL_3(QMDD, qmdd_gate, QMDD, gate_id_t, BDDVAR);

/* Applies given controlled gate to |q>. (Wrapper function) */
#define qmdd_cgate(qmdd,gate,c,t) (RUN(qmdd_cgate,qmdd,gate,c,t))
#define qmdd_cgate2(qmdd,gate,c1,c2,t) (RUN(qmdd_cgate2,qmdd,gate,c1,c2,t))
#define qmdd_cgate3(qmdd,gate,c1,c2,c3,t) (RUN(qmdd_cgate3,qmdd,gate,c1,c2,c3,t))
TASK_DECL_4(QMDD, qmdd_cgate,  QMDD, gate_id_t, BDDVAR, BDDVAR);
TASK_DECL_5(QMDD, qmdd_cgate2, QMDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR);
TASK_DECL_6(QMDD, qmdd_cgate3, QMDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR, BDDVAR);

/* Applies given controlled gate to |q>. (Wrapper function) */
#define qmdd_cgate_range(qmdd,gate,c_first,c_last,t) (RUN(qmdd_cgate_range,qmdd,gate,c_first,c_last,t))
TASK_DECL_5(QMDD, qmdd_cgate_range, QMDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR);

/**
 * Recursive implementation of applying single qubit gates
 */
#define qmdd_gate_rec(q,gate,target) (RUN(qmdd_gate_rec,q,gate,target))
TASK_DECL_3(QMDD, qmdd_gate_rec, QMDD, gate_id_t, BDDVAR);

/**
 * Recursive implementation of applying controlled gates
 */
#define qmdd_cgate_rec(q,gate,cs,t) (RUN(qmdd_cgate_rec,q,gate,cs,0,t))
TASK_DECL_5(QMDD, qmdd_cgate_rec, QMDD, gate_id_t, BDDVAR*, uint32_t, BDDVAR);

/**
 * Recursive implementation of applying controlled gates where the controlles 
 * are defined by a range 'c_first' through 'c_last'.
 */
#define qmdd_cgate_range_rec(q,gate,c_first,c_last,t) (RUN(qmdd_cgate_range_rec,q,gate,c_first,c_last,t,0))
TASK_DECL_6(QMDD, qmdd_cgate_range_rec, QMDD, gate_id_t, BDDVAR, BDDVAR, BDDVAR, BDDVAR);

/******************************</Applying gates>*******************************/





/*********************<Applying (controlled) sub-circuits>*********************/

// Circuit IDs
typedef enum circuit_id {
    CIRCID_swap,
    CIRCID_reverse_range,
    CIRCID_QFT,
    CIRCID_QFT_inv
} circuit_id_t;

/**
 * Circuit which implements a SWAP gate from single-qubit and controlled gates.
 */
QMDD qmdd_circuit_swap(QMDD qmdd, BDDVAR qubit1, BDDVAR qubit2);

/**
 * Circuit which reverses the order of the qubits in the given range.
 */
QMDD qmdd_circuit_reverse_range(QMDD qmdd, BDDVAR first, BDDVAR last);

/**
 * Executes the QFT circuit on qubits `first` through `last`.
 */
QMDD qmdd_circuit_QFT(QMDD qmdd, BDDVAR first, BDDVAR last);

/**
 * Executes the inverse QFT circuit on qubits `first` through `last`.
 */
QMDD qmdd_circuit_QFT_inv(QMDD qmdd, BDDVAR first, BDDVAR last);

/**
 * Applies the given circuit (parameters can be two targets or a range
 * depending on the circuit.)
 */
QMDD qmdd_circuit(QMDD qmdd, circuit_id_t circ_id, BDDVAR t1, BDDVAR t2);

/**
 * Generalized implementation of applying controlled versions of sub-circuit
 * functions defined here.
 * @param circ_id CIRCID_something
 * @param cs BDDVAR[] of control qubits. Needs to be length 3. If using fewer
 *           controls use e.g. cs = [c1, c2, AADD_INVALID_VAR]
 * @param t1 BDDVAR. Parameter 1 for given circuit.
 * @param t2 BDDVAR. Parameter 2 for given circuit.
 */
#define qmdd_ccircuit(qmdd, circ_id, cs, t1, t2) (RUN(qmdd_ccircuit,qmdd,circ_id,cs,0,t1,t2));
TASK_DECL_6(QMDD, qmdd_ccircuit, QMDD, circuit_id_t, BDDVAR*, uint32_t, BDDVAR, BDDVAR);

/**
 * Applies a phase of -1 to a single basis state |x>.
 * This is a CZ gate where we control on all qubits and when x_k = 0 we control
 * qubit k on 0, and where x_k = 1 we control qubit k on 1.
 * 
 * @param qmdd A QMDD encoding some quantum state |psi>.
 * @param n Number of qubits.
 * @param x A bitstring x of some computational basis state |x>.
 * 
 * TODO: generalize this to control on some but not all qubits.
 */
QMDD qmdd_all_control_phase(QMDD qmdd, BDDVAR n, bool *x);

/********************</Applying (controlled) sub-circuits>*********************/





/*******************************<Miscellaneous>********************************/

/**
 * Removes any global phase from state vector, assumes <qmdd> is the root edge.
 */
QMDD qmdd_remove_global_phase(QMDD qmdd);

/**
 * Trigger for gc of node table every n gates
 */
void qmdd_set_periodic_gc_nodetable(int every_n_gates);

/******************************</Miscellaneous>********************************/





/*******************************<Logging stats>********************************/

void qmdd_stats_start(FILE *out);
void qmdd_stats_set_granularity(uint32_t g); // log every 'g' gates (default 1)
void qmdd_stats_log(QMDD qmdd);
uint64_t qmdd_stats_get_nodes_peak();
double qmdd_stats_get_nodes_avg();
uint64_t qmdd_stats_get_logcounter();
void qmdd_stats_finish();

/******************************</Logging stats>********************************/





/****************************<Debug functionality>*****************************/

/**
 * Turns on/off (expensive) sanity checks
 */
void qmdd_set_testing_mode(bool on);

/**
 * Brute force sanity check to verify a given QMDD encodes (something close to) 
 * a unit vector.
 */
bool qmdd_is_close_to_unitvector(QMDD qmdd, BDDVAR n, double tol);
bool qmdd_is_unitvector(QMDD qmdd, BDDVAR n);

/**
 * Get the magnitude of a given qmdd. (Should equal 1 if the qmdd represents a
 * (pure) quantum state.
*/
double qmdd_get_magnitude(QMDD qmdd, BDDVAR n);

/***************************</Debug functionality>*****************************/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
