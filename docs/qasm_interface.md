# Q-Sylvan QASM interface

Q-Sylvan's QASM simulator currently supports a subset of [OpenQASM 2.0](https://en.wikipedia.org/wiki/OpenQASM). The supported functionality is outlined below.

##  Initialization
Currently, only one quantum and one classical register can be used.
* `qreg q[n];` : Declares an n-qubit quantum register `q`.
* `creg c[n];` : Declares an n-bit classical register `c`.


## Gate operations
For single qubit gates on qubit k:
* `i q[k];` : Identity gate
* `x q[k];` : Pauli-X gate
* `y q[k];` : Pauli-Y gate
* `z q[k];` : Pauli-Z gate
* `s q[k];` : S gate (= &radic;Z)
* `t q[k];` : T gate (= &radic;S)
* `h q[k];` : Hadamard gate
* `sx q[k];` : &radic;X gate
* `sy q[k];` : &radic;Y gate
* `sdg q[k];` : S<sup>&dagger;</sup> gate
* `tdg q[k];` : T<sup>&dagger;</sup> gate
* `rx(a) q[k];` : Rotation around x-axis with angle 2&pi; &times; a
* `ry(a) q[k];` : Rotation around y-axis with angle 2&pi; &times; a
* `rz(a) q[k];` : Rotation around z-axis with angle 2&pi; &times; a

For controlled gates we can prefix a `c` to an existing gate, e.g. for a CX gate:
* `cx q[c], q[t];` : Controlled X gate on control c and target t<sup>[1](#control_order)</sup>.

This can also also be extended to arbitrarily many control qubits, where the last argument is the target qubit, and the others the controls:
* `cccx q[c1], q[c2], q[c3], q[t];` Controlled X gate with 3 control qubits<sup>[1](#control_order)</sup>.

<sub><a id="control_order">1</a>: These control and targets can be in any order, but the overhead is the least when all c<sub>k</sub> < t.</sub>


## Measurements
* `measure q[k]->c[j];` : Computational basis measurement of qubit k. The measurement outcome is stored in `c[j]`.


## Classical conditioning
There are two ways to perform a quantum gate conditioned on a classical expression.
* `if(c[k]) x q[0];` : perform an X gate on qubit 0 if `c[k]` holds the value 1.
* `if(c==5) x q[0];` : perform an X gate on qubit 0 if the bits in the classical register represent the number 5


## Other
* `// comment` : Any line starting with `//` is ignored
* `barrier;` : Prevents transformations of the circuit across this column. When the simulator is set to "balanced" mode barriers also switch the the computation from matrix-matrix to matrix-vector and vice versa.
