import matplotlib.pyplot as plt
from numpy.core.records import array
from numpy.testing._private.utils import assert_equal
import pandas as pd
from math import pi
import numpy as np
import subprocess
import itertools
import random
import math
import time
import sys

import cirq

# Initialisation of some parameters
seed = 8092             # Random seed
np.random.seed(seed)    # Initialising andom seed

# Builds the variational circuit
def circuit(q, W, theta, layers):
    for _ in range(2):
        for i in range(len(q)):
            yield cirq.H(q[i])
        for i in range(len(q)):
            yield cirq.ZPowGate(exponent=W[i]/pi)(q[i])
        for i in range(len(q)-1):
            for j in range(i+1,len(q)):
                yield cirq.ZPowGate(exponent=((pi-W[i])*(pi-W[j]))/pi).on(q[j]).controlled_by(q[i])
    for i in range(len(q)):
        yield cirq.ZPowGate(exponent = theta[2*i]/pi)(q[i])
        yield cirq.ry(theta[2*i + 1])(q[i])
    for j in range(1,layers+1):
        t = theta[range(6*(j),6*(j+1))]
        for i in range(len(q)):
            yield cirq.CZ.on(q[(i+1)%len(q)],q[i])
        for i in range(len(q)):
            yield cirq.ZPowGate(exponent=t[2*i]/pi)(q[i])
            yield cirq.ry(t[2*i+1])(q[i])
    for i in range(len(q)):
        yield cirq.measure(q[i], key=str(i))

class Circuit():
    def __init__(self, qubits, layers, gates):
        self.qubits = qubits
        self.layers = layers
        self.gates = gates
        self.circuit_file = "./vqc_qasm.txt"
        self.unresolved_circuit = self.create_circuit()
        self.resolved_circuit = None

    def create_circuit(self):
        # Initial QASM setup
        QASM_code = ['OPENQASM 2.0;']
        QASM_code.append('include "qelib1.inc";')
        QASM_code.append('qreg q[{}];'.format(self.qubits))
        QASM_code.append('creg c[{}];'.format(self.qubits))
        # Entangle datapoint
        for _ in range(2):
            for i in range(self.qubits):
                QASM_code.append('h q[{}];'.format(i))
            for i in range(self.qubits):
                QASM_code.append('rz({})'+' q[{}];'.format(i))
            for (i,j) in itertools.combinations(range(self.qubits), 2):
                QASM_code.append('crz({})'+' q[{}], q[{}];'.format(i, j))
        # Parametrized gates
        for gate in self.gates:
            for i in range(qubits):
                QASM_code.append('{}'.format(gate)+'({})'+' q[{}];'.format(i))
        for _ in range(self.layers):
            for i in range(self.qubits-1):
                QASM_code.append('cz q[{}], q[{}];'.format(i, i+1))
            if self.qubits != 2:
                QASM_code.append('cz q[{}], q[{}];'.format(0, self.qubits-1))
            for gate in self.gates:
                for i in range(qubits):
                    QASM_code.append('{}'.format(gate)+'({})'+' q[{}];'.format(i))
        for i in range(qubits):
            QASM_code.append('measure q[{}]->c[{}];'.format(i, i))
        return QASM_code

    def resolve_parameters(self, X, theta):
        self.resolved_circuit = self.unresolved_circuit.copy()
        j = 0
        # Create entangled rotations for every qubit pair
        X_ent = [(1-X[item[0]])*(1-X[item[1]]) for item in itertools.combinations(range(self.qubits),2)]
        params = np.concatenate((X,X_ent,X,X_ent,theta))
        for i in range(len(self.resolved_circuit)):
            if(self.resolved_circuit[i].find('{}') != -1):
                self.resolved_circuit[i] = self.resolved_circuit[i].format(round(params[j]%1,3))
                j += 1

    def run(self, X, theta, shots):
        # Resolve parameters
        self.resolve_parameters(X, theta)
        f = open(self.circuit_file, "w")
        for line in self.resolved_circuit:
                f.write(line+'\n')
        f.close()
        # Run circuit
        output = subprocess.run(["../build/parser/QASM_to_Sylvan", "-f", self.circuit_file, "-s", str(shots),"-r",str(int(random.randint(0,sys.maxsize)))], stdout=subprocess.PIPE)
        output = output.stdout.decode('utf-8')
        return output.split("\n")[:-1]


class VQC:
    def __init__(self, qubits, layers, gates=["ry"], shots=1000):
        self.qubits = qubits
        self.layers = layers
        self.gates = gates
        self.shots = shots
        self.nr_params = (self.qubits*len(self.gates))*(self.layers+1)
        self.theta = np.random.rand(self.nr_params,)
        self.theta = np.append(self.theta,0)
        self.a = 2.5
        self.c = 0.1
        self.z = self.a/24
        self.alpha = 0.602
        self.gamma = 0.101
        self.circuit = Circuit(self.qubits, self.layers, self.gates)

    def train(self, X, Y, iterations):
        start = end = time.time()
        total_loss = np.zeros(iterations)
        # Start iterations
        for k in range(1,iterations+1):
            print('Finished iteration',k,'/',iterations,'at: ',end='')
            # Update parameters
            c_n = self.c/(k**(self.gamma))
            a_n = self.a/(k**(self.alpha))
            z_n = self.z/(k**(self.alpha))
            delta_n = 2*np.random.randint(2, size = self.nr_params+1) - 1
            # Run variational classifier with theta+delta and theta-delta
            p_plus  = self.J(self.theta+c_n*delta_n, X)
            p_minus = self.J(self.theta-c_n*delta_n, X)
            # Calculate the loss for each run
            loss_plus  = self.empirical_risk(Y, p_plus) 
            loss_minus = self.empirical_risk(Y, p_minus)

            # Compute gradient and update theta accordingly
            grad = delta_n*((loss_plus - loss_minus)/(2*c_n))
            self.theta[:-1] = (self.theta[:-1] - a_n*grad[:-1]) #% (2*pi)
            self.theta[-1] = (self.theta[-1] - z_n*grad[-1]) 
            # parameter b is probably taking too large steps.

            # Save average loss for plotting
            total_loss[k-1] = (loss_plus + loss_minus)/2
            end = time.time()
            print(np.round(end-start,3),'seconds')
            print('loss: ', total_loss[k-1])
        # Print time taken for all iterations
        return total_loss

    def test(self, X, Y):
        p_t= self.J(self.theta, X)
        acc_t = self.accuracy(Y, p_t, self.theta[-1])
        print("accuracy: ", acc_t)

    # Runs the circuit and parses the result to a probability estimate
    def J(self, gradient_theta, X):
        p = np.zeros([len(X),])
        for i in range(len(X)):
            results = self.circuit.run(X[i], gradient_theta, shots=self.shots)
            p[i] = self.probability_estimate(results)
        return p

    # Runs the circuit and parses the result to a probability estimate
    def old_J(self, gradient_theta, X):
        q = [cirq.GridQubit(i, 0) for i in range(self.qubits)]
        simulator = cirq.Simulator()
        p = np.zeros([len(X),])
        for i in range(len(X)):
            c = cirq.Circuit()
            c.append(circuit(q, X[i], gradient_theta, self.layers))
            results = simulator.run(c, repetitions=self.shots)
            p[i] = self.old_probability_estimate(results)
        return p

    # Returns a probability of seeing label 1 for a datapoint
    def probability_estimate(self, results):
        final_prob = 0
        for var in results:
            measurement, prob = var.split(": ")
            if measurement.count('1') % 2 == 1:
                final_prob += int(prob)
        return final_prob/self.shots

    # Returns a probability of seeing label 1 for a datapoint
    def old_probability_estimate(self, results):
        counter = (results.multi_measurement_histogram(keys="012"))
        p_hold = 0
        for j in counter:
            if j.count(1) % 2 == 1:
                p_hold += counter[j] 
        return p_hold/self.shots

    # Returns the emperical risk of the predictions made by the classifier
    def empirical_risk(self, labels, p):
        R = 200
        x = math.sqrt(R)*(np.multiply(labels,0.5-p) + 0.5*labels*self.theta[-1])/(np.sqrt(2*np.multiply(p,1-p))+1e-10)
        loss = (1 / (1 + np.exp(-x)))
        return sum(loss) / len(labels)

    # Returns the accuracy over the predicted labels of a dataset
    def accuracy(self, labels, p_hat, b):
        est_label = self.assign_label(p_hat, b)
        err = (labels*est_label - 1)/-2
        acc = 1 - np.sum(err)/len(err)
        return acc

    # Assigns labels based on prediction and bias
    def assign_label(self, p, b):
        labels = np.ones([len(p),])*-1
        for i in range(len(p)):
            if (p[i] > ((1 - p[i]) - b)):
                labels[i] = 1
        return labels

def get_data(filename, variables):
    # Load the data and split parameters and labels
    df = pd.read_csv(filename)
    X = df.iloc[:,:variables].to_numpy()
    Y = df.iloc[:,variables].to_numpy()
    Y = 2*Y - 1
    # Initialise training data (make sure it is balanced by repeating until it is)
    get_data = True
    ratio = sum((Y+1)/2)/len(Y)
    while get_data:
        random.shuffle(X)
        random.shuffle(Y)
        X_train = X[:int(np.ceil(len(X)*4/5))]
        X_test = X[int(np.floor(len(X)*4/5)):]
        Y_train = Y[:int(np.ceil(len(Y)*4/5))]
        Y_test = Y[int(np.floor(len(Y)*4/5)):]
        ratio_t = sum((Y_train+1)/2)/len(Y_train)
        if ratio-0.02 <= ratio_t and ratio_t <= ratio+0.02:
            return (X_train, Y_train), (X_test, Y_test)

if __name__ == "__main__":
    qubits = 3
    layers = 4
    gates = ["ry", "rz"]
    shots = 2000
    iterations = 50

    vqc = VQC(qubits, layers, gates, shots)
    data_train, data_test = get_data("QA_data_2pi_unbal.csv", variables=qubits)
    loss = vqc.train(X=data_train[0], Y=data_train[1], iterations=iterations)
    print(loss)
    vqc.test(X=data_train[0], Y=data_train[1])
    vqc.test(X=data_test[0], Y=data_test[1])
