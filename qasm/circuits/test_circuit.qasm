OPENQASM 2.0;
include "qelib1.inc";

// 15 qubit quantum register and 2 bit classical register
qreg q[15];
creg c[15];

// Create |Phi^+>
x q[0];
y q[1];
z q[3];
h q[4];
s q[5];
sdg q[6];
t q[7];
tdg q[8];
sx q[9];
sxdg q[10];
rx(pi) q[11];
ry(pi) q[12];
rz(pi) q[13];
p(pi) q[14];


// Measure state

