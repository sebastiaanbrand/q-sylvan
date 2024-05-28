"""
Testing running sim_qasm from command line.
"""
import os
import subprocess
import json
import numpy as np


TOLERANCE = 1e-6
SIM_QASM = './build/qasm/sim_qasm'
QASM_DIR = 'qasm/test/test_circuits/'


def fidelity(a, b):
    """
    Compute the fidelity of two vectors.
    """
    in_prod = np.dot(a.conj().T, b)
    fidelity = (abs(in_prod))**2
    return fidelity


def get_vector(qasm_file : str):
    """
    Simulate given quantum circuit and return state vector.
    """
    output = subprocess.run([SIM_QASM, '--state-vector', qasm_file], 
                            stdout=subprocess.PIPE)
    data = json.loads(output.stdout)
    vector = np.apply_along_axis(lambda args: [complex(*args)], 1,
                                 data['state_vector']).flatten()
    return vector


def test_bell_psi_plus():
    """
    Test |Psi+> = 1/sqrt(2)[1 0 0 1]
    """
    qasm_file = os.path.join(QASM_DIR, 'bell_psi_plus.qasm')
    vector = get_vector(qasm_file)
    reference = np.array([1/np.sqrt(2) + 0j, 
                          0 + 0j,
                          0 + 0j,
                          1/np.sqrt(2) + 0j])
    assert abs(fidelity(vector, reference) - 1) < TOLERANCE


def test_bell_psi_min():
    """
    Test |Psi-> = 1/sqrt(2)[1 0 0 -1]
    """
    qasm_file = os.path.join(QASM_DIR, 'bell_psi_min.qasm')
    vector = get_vector(qasm_file)
    reference = np.array([1/np.sqrt(2) + 0j, 
                          0 + 0j,
                          0 + 0j,
                          -1/np.sqrt(2) + 0j])
    assert abs(fidelity(vector, reference) - 1) < TOLERANCE


def test_bell_phi_plus():
    """
    Test |Phi+> = 1/sqrt(2)[0 1 1 0]
    """
    qasm_file = os.path.join(QASM_DIR, 'bell_phi_plus.qasm')
    vector = get_vector(qasm_file)
    reference = np.array([0 + 0j, 
                          1/np.sqrt(2) + 0j,
                          1/np.sqrt(2) + 0j,
                          0 + 0j])
    assert abs(fidelity(vector, reference) - 1) < TOLERANCE


def test_bell_phi_min():
    """
    Test |Phi-> = 1/sqrt(2)[0 1 -1 0]
    """
    qasm_file = os.path.join(QASM_DIR, 'bell_phi_min.qasm')
    vector = get_vector(qasm_file)
    reference = np.array([0 + 0j, 
                          1/np.sqrt(2) + 0j,
                          -1/np.sqrt(2) + 0j,
                          0 + 0j])
    assert abs(fidelity(vector, reference) - 1) < TOLERANCE


