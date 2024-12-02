from qiskit import QuantumCircuit
from qiskit.qasm2 import dumps
import random
import os

import globals_pipeline as g

# Generates a random quantum circuit based on Clifford gates and T-gate.
def generate_random_clifford_t_circuit(num_qubits, num_gates):
 
    # Available gates
    single_qubit_clifford_gates = ['h', 's', 'sdg', 'x', 'y', 'z']  
    two_qubit_clifford_gates = ['cx', 'cz']
    t_gate = ['t'] 

    # Initialize quantum circuit
    qc = QuantumCircuit(num_qubits)

    single_count = 0
    two_count = 0
    t_count = 0

    for _ in range(num_gates):

        gate_types = ['single', 'two', 't']

        probability_gate_types = [0.5 * (1 - g.PERCENTAGE_T_GATES / 100), 0.5 * (1 - g.PERCENTAGE_T_GATES / 100), g.PERCENTAGE_T_GATES / 100]

        gate_type = random.choices(gate_types, weights=probability_gate_types, k=1)
        
        if gate_type[0] == 'single':

            single_count = single_count + 1

            # Apply a single-qubit Clifford gate to a random qubit
            qubit = random.randint(0, num_qubits - 1)
            gate = random.choice(single_qubit_clifford_gates)

            if gate == 'h':
                qc.h(qubit)
            elif gate == 's':
                qc.s(qubit)
            elif gate == 'sdg':
                qc.sdg(qubit)
            elif gate == 'x':
                qc.x(qubit)
            elif gate == 'y':
                qc.y(qubit)
            elif gate == 'z':
                qc.z(qubit)

        elif gate_type[0] == 'two' and num_qubits > 1:

            two_count = two_count + 1

            # Apply a two-qubit Clifford gate to a random pair of qubits
            control = random.randint(0, num_qubits - 1)
            target = random.randint(0, num_qubits - 1)

            # Ensure target and control are different
            while target == control:  
                target = random.randint(0, num_qubits - 1)

            gate = random.choice(two_qubit_clifford_gates)
            if gate == 'cx':
                qc.cx(control, target)
            elif gate == 'cz':
                qc.cz(control, target)

        elif gate_type[0] == 't':

            t_count = t_count + 1

            # Apply a T-gate to a random qubit
            qubit = random.randint(0, num_qubits - 1)
            qc.t(qubit)

    #print(f"Number of single gates = {single_count}")
    #print(f"Number of two gates    = {two_count}")
    #print(f"Number of t gates      = {t_count}")
    #print(f"Total gates            = {single_count + two_count + t_count}")

    return qc

def save_circuit_as_qasm(circuit, file_name):
    with open(file_name, 'w') as file:
        file.write(dumps(circuit))

# Function to generate and save multiple circuits
def generate_and_save_circuits(num_qubits, num_gates):

    # Ensure the output directory exists
    os.makedirs(g.QASM_DIR, exist_ok=True)

    for i in range(g.NUM_CIRCUITS):
        circuit = generate_random_clifford_t_circuit(num_qubits, num_gates)
        file_name = os.path.join(g.QASM_DIR, f"clifford_T_circuit_{i+1}_{num_qubits}_{num_gates}.qasm")
        save_circuit_as_qasm(circuit, file_name)
        print(f"Circuit {i+1} saved to {file_name}")



