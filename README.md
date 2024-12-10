# Q-Sylvan

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![CI testing](https://github.com/sebastiaanbrand/q-sylvan/actions/workflows/cmake.yml/badge.svg)](https://github.com/sebastiaanbrand/q-sylvan/actions/workflows/cmake.yml)

Q-Sylvan extends the parallel decision diagram library [Sylvan](https://github.com/trolando/sylvan) (v1.8.0) with QMDDs (i.e. factored EVBDDs with complex edge weights), as well as functionality to simulate quantum circuits. This is currently still a beta version.


## Installation

### Dependencies

Q-Sylvan requires the following libraries: `popt`, `GMP`, `MPFR` and `MPC`. On Ubuntu it should be possible to install these with
- `sudo apt-get install libpopt-dev`
- `sudo apt-get install libgmp-dev`
- `sudo apt-get install libmpfr-dev`
- `sudo apt-get install libmpc-dev`
- `sudo apt-get install cmake`


### Compiling the code
After downloading or cloning the repository, from the repository folder the code can be compiled with:
1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

Test can be run with `ctest`.

## Simulations of circuits presented in qasm on QMDD and MTBDD

Running can be activated with `./run_qasm_on_qmdd <filename>.qasm --state-vector` in `build/qasm/` directory.

For MTBDD type `./run_qasm_on_mtbdd <filename>.qasm --state-vector` in `build/qasm` directory.

Debugging can be activated with `gdb ./run_qasm_on_qmdd` in `build/qasm` directory.

The GNU GDB Debugger accepts commands as `(gdb) run <filename>.qasm --state-vector`.

## Tests of the simulators based on QMDD and MTBDD

Running tests can be activated with `pytest -v` in the root directory, `Q-Sylvan/q-sylvan/`.


## Example usage
The following code snippets are a toy example for creating the [Bell state](https://en.wikipedia.org/wiki/Bell_state) `|Phi^+> = 1/sqrt(2) (|01> + |10>)`, first in the C interface, and secondly in the QASM interface.

### C interface
```C
// Create |Phi^+>
int nqubits = 2;
QMDD state = qmdd_create_all_zero_state(nqubits);
state = qmdd_gate(state, GATEID_H, 0);     // H on q0
state = qmdd_cgate(state, GATEID_X, 0, 1); // CNOT on c=q0, t=q1
state = qmdd_gate(state, GATEID_X, 0);     // X on q0

// Measure state
bool outcome[] = {0, 0};
double prob;
qmdd_measure_all(state, nqubits, outcome, &prob);
```
This code can be found in [`examples/bell_state.c`](examples/bell_state.c) and after compiling the code as described above can be run with `./examples/bell_state` from the `build/` directory. A more complete set of available functions can be found [here](docs/documentation/c_interface.md).

### QASM interface
```C
OPENQASM 2.0;
include "qelib1.inc";

// 2 qubit quantum register and 2 bit classical register
qreg q[2];
creg c[2];

// Create |Phi^+>
h q[0];
cx q[0], q[1];
x q[0];

// Measure state
measure q[0]->c[0];
measure q[1]->c[1];
```
This code can be found in [`qasm/circuits/bell_state.qasm`](qasm/circuits/bell_state.qasm) and can be run with `./build/qasm/run_qasm_on_qmdd qasm/circuits/bell_state.qasm`. All gates specified in [qelib1.inc](https://github.com/Qiskit/qiskit-terra/blob/main/qiskit/qasm/libs/qelib1.inc) are supported, as well as measurements. Custom definitions of gates and classical conditioning are currently not supported.


## Documentation
A more complete documentation of the C interface can be found [here](docs/documentation/c_interface.md).


## Acknowledgements
This work is supported by the [NEASQC](https://cordis.europa.eu/project/id/951821) project, funded by the European Union's Horizon 2020 programme, Grant Agreement No. 951821.
