# Q-Sylvan

Q-Sylvan extends [Sylvan](https://github.com/trolando/sylvan) with QMDDs, as well as functionality to simulate quantum circuits.


## Installation

### Dependencies
Q-Sylvan requires the following libraries: `popt` and `GMP`. On Ubuntu it should be possible to install these with
- `sudo apt-get install libpopt-dev`
- `sudo apt-get install libgmp3-dev`


### Compiling the code
After downloading or cloning the repository, from the repository folder the code can be compiled with:
1. `mkdir build`
2. `cd build`
3. `cmake ..`
4. `make`

Test can be run with `make test`.

Installing should also be possible with `make install` but this will clash with any other version of Sylvan already installed.


## Example usage
The following code snippets are a toy example for creating the [Bell state](https://en.wikipedia.org/wiki/Bell_state) `|Phi^+> = 1/sqrt(2) (|01> + |10>)`, first in the C interface, and secondly in the QASM interface.

### C interface
```C
// Create |Phi^+>
int nqubits = 2;
QDD state = qdd_create_all_zero_state(nqubits);
state = qdd_gate(state, GATEID_H, 0);     // H on q0
state = qdd_cgate(state, GATEID_X, 0, 1); // CNOT on c=q0, t=q1
state = qdd_gate(state, GATEID_X, 0);     // X on q0

// Measure state
bool outcome[] = {0, 0};
double prob;
qdd_measure_all(state, nqubits, outcome, &prob);
```
This code can be found in `examples/bell_state.c` and after compiling the code as described above can be run with `./examples/bell_state` from the `build/` directory.

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
This code can be found in `qasm/circuits/bell_state.qasm` and can be run with `./qasm/circuit_to_Sylvan ../qasm/circuits/bell_state.qasm -r=100` from the `build/` directory.


## Folders overview
```
TODO
```

