"""
Compute state vector with qiskit. Used for generating tests.
"""
import argparse
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

parser = argparse.ArgumentParser(description='Compute reference state using qiskit.')
parser.add_argument('qasm_file', help='Path to .qasm file.')


def main():
    args = parser.parse_args()
    qc = QuantumCircuit.from_qasm_file(args.qasm_file)
    qc.remove_final_measurements()
    print(qc)
    statevector = Statevector(qc)
    print(repr(statevector.data))


if __name__ == '__main__':
    main()
