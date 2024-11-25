from qiskit import QuantumCircuit
from qiskit.qasm2 import dumps
import random
import os

# Generates a random quantum circuit based on Clifford gates and T-gate.
def generate_random_clifford_t_circuit(num_qubits, num_gates):
 
    # Available gates
    single_qubit_clifford_gates = ['h', 's', 'sdg', 'x', 'y', 'z']  
    two_qubit_clifford_gates = ['cx', 'cz', 'swap', 'iswap'] 
    t_gate = ['t'] 

    # Initialize quantum circuit
    qc = QuantumCircuit(num_qubits)

    for _ in range(num_gates):
        gate_type = random.choice(['single', 'two', 't'])
        
        if gate_type == 'single':
            
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

        elif gate_type == 'two' and num_qubits > 1:
            
            # Apply a two-qubit Clifford gate to a random pair of qubits
            control = random.randint(0, num_qubits - 1)
            target = random.randint(0, num_qubits - 1)
            while target == control:  # Ensure target and control are different
                target = random.randint(0, num_qubits - 1)
            gate = random.choice(two_qubit_clifford_gates)
            if gate == 'cx':
                qc.cx(control, target)
            elif gate == 'cz':
                qc.cz(control, target)
            elif gate == 'swap':
                qc.swap(control, target)
            elif gate == 'iswap':
                qc.iswap(control, target)
        
        elif gate_type == 't':

            # Apply a T-gate to a random qubit
            qubit = random.randint(0, num_qubits - 1)
            qc.t(qubit)

    return qc

def save_circuit_as_qasm(circuit, file_name):
    with open(file_name, 'w') as file:
        file.write(dumps(circuit))

# Function to generate and save multiple circuits
def generate_and_save_circuits(num_circuits, num_qubits, num_gates, output_dir):

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for i in range(num_circuits):
        circuit = generate_random_clifford_t_circuit(num_qubits, num_gates)
        file_name = os.path.join(output_dir, f"clifford_T_circuit_{i+1}.qasm")
        save_circuit_as_qasm(circuit, file_name)
        print(f"Circuit {i+1} saved to {file_name}")

# Parameters
num_circuits = 50   # Number of circuits to generate
num_qubits = 100    # Number of qubits per circuit
num_gates = 1000    # Number of gates per circuit
output_dir = "qasm_clifford_T"  # Directory to save the circuits

# Generate and save circuits
generate_and_save_circuits(num_circuits, num_qubits, num_gates, output_dir)
