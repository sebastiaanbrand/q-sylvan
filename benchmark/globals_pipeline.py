#
# Globals
#

global NUM_CIRCUITS
global MIN_NUM_QUBITS
global MAX_NUM_QUBITS
global STEPSIZE_NUM_QUBITS
global MIN_NUM_GATES
global MAX_NUM_GATES
global STEPSIZE_NUM_GATES
global PERCENTAGE_T_GATES

global PRECISION
global METHOD
global SIM_QASM_EXE

global QASM_DIR
global JSON_DIR
global PLOTS_DIR

#
# Parameters
#

NUM_CIRCUITS = 1             # Number of random circuits to generate per number of qubits and number of gates
MIN_NUM_QUBITS = 10           # Min number of qubits per circuit
MAX_NUM_QUBITS = 200          # Max number of qubits per circuit
STEPSIZE_NUM_QUBITS = 10     # Stepsize of the number of qubits
MIN_NUM_GATES = 100          # Min number of gates per circuit
MAX_NUM_GATES = 1000         # Max number of gates per circuit
STEPSIZE_NUM_GATES = 100     # Stepsize of the number of gates
PERCENTAGE_T_GATES = 55      # Percentage T gates in quantum circuit

#PRECISION = 16
#METHOD = 'MTBDD'
#SIM_QASM_EXE = './../build/qasm/run_qasm_on_mtbdd'

PRECISION = 64
METHOD = 'QMDDmax'
SIM_QASM_EXE = './../build/qasm/run_qasm_on_qmdd'

QASM_DIR  = f'qasm_clifford_T_mixed_{METHOD}_{PERCENTAGE_T_GATES}_percent_{MIN_NUM_QUBITS}_{MAX_NUM_QUBITS}_{MIN_NUM_GATES}_{MAX_NUM_GATES}/'
JSON_DIR  = f'json_clifford_T_mixed_{METHOD}_{PERCENTAGE_T_GATES}_percent_{MIN_NUM_QUBITS}_{MAX_NUM_QUBITS}_{MIN_NUM_GATES}_{MAX_NUM_GATES}/'
PLOTS_DIR = f'plot_clifford_T_mixed_{METHOD}_{PERCENTAGE_T_GATES}_percent_{MIN_NUM_QUBITS}_{MAX_NUM_QUBITS}_{MIN_NUM_GATES}_{MAX_NUM_GATES}/'

