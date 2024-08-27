# Q-Sylvan C interface functions

NOTE: In the code we are refering to the QMDDs as QMDDs, but our QMDDs are really just QMDDs.

## Initialization of states vectors / matrices
* `qmdd_create_all_zero_state(int n)` : Creates a QMDD for an n-qubit state |00...0>. 
* `qmdd_create_basis_state(int n, bool* x)` : Creates a QMDD for an n-qubit state |x>.
* `qmdd_create_all_identity_matrix(int n)` : Creates a QMDD representing an n-qubit identity matrix.
* `qmdd_create_single_qubit_gate(int n, int t, gate_id_t gateid)` : Creates an n-qubit matrix QMDD which applies gate `gateid` to qubit t and I to all others.
* `qmdd_create_single_qubit_gates(int n, gate_id_t *gateid)` : Creates an n-qubit matrix QMDD which applies the given list of n `gatesid`'s to n qubits.
* `qmdd_create_single_qubit_gates_same(int n, gate_id_t gateid)` : Creates an n-qubit matrix QMDD which applies single qubit gate `gateid` to all qubits.
* `qmdd_create_cgate(int n, int c, int t, gate_id_t gateid)` : Creates an n-qubit matrix QMDD which acts as a controlled-`gateid` gate on target qubit t with control qubit c, and I on all others. (Requires c < t.)
* `qmdd_create_multi_cgate(int n, int *c_options, gate_id_t gateid)` : Creates a controlled-`gateid` gate which acts on the qubits as specified in by the `c_options` array. `c_options[k]` can contain `-1` = ignore qubit k (apply I), `0` = control qubit k on |0>, `1` = control qubit k on |1>, and `2` = the target qubit. (Requires the index of the target qubit to be greater than the indices of the controls.)
* `qmdd_create_all_control_phase(int n, bool *x)` : Creates an n-qubit CZ gate, controlled on all qubits. The array x specifies whether each qubit k is controlled on |0> (when x[k]=0) or on |1> (when x[k]=1).

## Gate operations
* `qmdd_gate(QMDD qmdd, gate_id_t gateid , int t)` : Applies given (single qubit) gate to qubit t.
* `qmdd_cgate(QMDD qmdd, gate_id_t gateid, int c, int t)` : Applies controlled-`gateid` gate to (c)ontrol and (t)arget. (Requires c < t.)
* `qmdd_cgate2(QMDD qmdd, gate_id_t gateid, int c1, int c2, int t)` : As above but with two controls (c1 < c2 < t).
* `qmdd_cgate3(QMDD qmdd, gate_id_t gateid, int c1, int c2, int c3, int t)` : As above but with three controls (c1 < c2 < c3 < t).
* `qmdd_cgate_range(QMDD qmdd, gate_id_t gateid , int c_first, int c_last, int t)` : Applies controlled-`gateid` to (t)arget, with all qubits between (and including) c_first and c_last as controls (c_first < c_last < t).
* `evbdd_matvec_mult(QMDD mat, QMDD vec, int n)` : Computes mat|vec> for an 2^n vector and a 2^n x 2^n matrix.
* `evbdd_matmat_mult(QMDD a, QMDD b, int)` : Computes a*b for two 2^n x 2^n matrices.
* `evbdd_vec_tensor_prod(QMDD a, QMDD b, int nqubits_a)` : Computes a \tensor b for two vector QMDDs.
* `evbdd_mat_tensor_prod(QMDD a, QMDD b, QMDD nqubits_a)` : Computes a \tensor b for two matrix QMDDs.

## Measurements and related
Note: For measurements, the post-measurement state is the return value of the measurement function, while the input qmdd is unaffected.
* `qmdd_get_amplitude(QMDD qmdd, bool* x, int nqubits)` : Get the amplitude <qmdd|x> as a complex struct, where x is a bool array of lenght n.
* `qmdd_measure_qubit(QMDD qqd, int k, int n, int *m, double *p)` : Measures qubit k of the given n-qubit state. `&m` will contain the measurement outcome, and `&p` will contain Pr(m = 0).
* `qmdd_measure_all(QMDD qmdd, int n, bool *ms, double *p)` : Does a (computational basis) measurement of all n qubits in the qmdd. The measurement outcomes are put in `ms`, which needs to be an bool array of lenght n. `p` (a pointer to a single double) will contain the probability |<qmdd|ms>|^2.

## Predefined gates
* `GATEID_I` : identity gate
* `GATEID_X` : Pauli-X gate
* `GATEID_Y` : Pauli-Y gate
* `GATEID_Z` : Pauli-Z gate
* `GATEID_H` : Hadamard gate
* `GATEID_S` : S gate (= &radic;Z)
* `GATEID_Sdag` : S<sup>&dagger;</sup> gate
* `GATEID_T` : T gate (= &radic;S)
* `GATEID_Tdag` : T<sup>&dagger;</sup> gate
* `GATEID_sqrtX` : &radic;X gate
* `GATEID_sqrtXdag` : (&radic;X)<sup>&dagger;</sup> gate
* `GATEID_sqrtY` : &radic;Y gate
* `GATEID_sqrtYdag` : (&radic;Y)<sup>&dagger;</sup> gate
* `GATEID_Rk(k)` : Phase gate with angle 2&pi; &times; 1/2<sup>k</sup>
* `GATEID_Rk_dag(k)` : Phase gate with angle 2&pi; &times; 1/2<sup>-k</sup>
* `GATEID_Rx(theta)` : Rotation around x-axis with angle theta
* `GATEID_Ry(theta)` : Rotation around y-axis with angle theta
* `GATEID_Rz(theta)` : Rotation around z-axis with angle theta


## Other
* `evbdd_countnodes(QMDD qmdd)` Counts the number of nodes in the given QMDD.
* `evbdd_fprintdot(FILE *out, QMDD qmdd, bool draw_zeroes)` Writes a .dot representation of the given QMDD to the given file.
