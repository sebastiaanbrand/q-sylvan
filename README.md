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

### C interface
The following code snippet creates the [Bell state](https://en.wikipedia.org/wiki/Bell_state) $\ket{\Phi^+} = \frac{1}{\sqrt{2}} (\ket{0} + \ket{1})$. This code can also be found in `examples/bell_state.c` and after compiling the code as described above can be run with `./examples/bell_state`.
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

### QASM interface


## Folders overview
```
TODO
```

