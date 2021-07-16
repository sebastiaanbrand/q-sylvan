#####################################################################################
#                                  How to use                                       #
# Make sure you fprintf to a .txt file in the folder "plotting"                     #
# The gates should be written to the file in one of two ways:                       #
#  1. A single gate should be represented as "q: G", where q is the qubit number    #
#        and G is the gate visual (using only one letter in capital case)           #
#  2. A controlled-gate should be represented as "c1/.../c2->t: G", where c1/.../c2 #
#        are all control-qubits, seperated by "/", t is the target qubitand G is    #
#        the gate visual (using only one letter in capital case)                    #
# Running this file, giving the .txt file as an argument, will print the circuit    #
# described in the given file                                                       #
#####################################################################################
import sys

if len(sys.argv) != 2:
    print("Pleas give exactly one filename.")
    exit(-1)
f = open(sys.argv[1], "r")

# Prepare text outputs
line = f.readline()
n = int(line.split(" ")[0])
text = ["---" for _ in range(int(n))]

barrier = False
# Parse gates
for line in f.readlines():
    line = line[:-1] # remove \n
    if line == "barrier":
        barrier = True
        for i in range(n):
            text[i] += ":-"
        continue
    qubits, command = line.split(": ")
    qubit_split = qubits.split("->")
    if(len(qubit_split) == 1): # single qubit gate
        q = int(qubit_split[0])
        if text[q][-2] == "-" and not barrier:
            t = text[q][:-2]
            text[q ]= t + command + "-"
        else:
            for i in range(n):
                if i == q: text[i] += command + "-"
                else: text[i] += "--"
    else: # controlled gate
        target = int(qubit_split[1])
        controls = qubit_split[0].split("/")
        controls = [int(i) for i in controls] 
        for i in range(n):
            if i in controls: text[i] += "@-"
            elif i == target:
                if command == "Z": text[i] += "@-"
                else: text[i] += command + "-"
            elif i > min(min(controls),target) and i < max(max(controls), target): text[i] += "|-"
            else: text[i] += "--"
    barrier = False
for line in text:
    print(line)