OPENQASM 2.0;
include "qelib1.inc";

// 16 qubit quantum register and 2 bit classical register
qreg q[4];
//creg c[4];

// Unitary gates
x q[0];
//y q[1];
//z q[3];
//h q[4];
//s q[5];
//sdg q[6];
//t q[7];
//tdg q[8];
//sx q[9];
//sxdg q[10];
//rx(pi) q[11];
//rx(2*pi) q[11];
//ry(pi) q[12];
//rz(pi) q[13];
//p(pi) q[14];
//u(pi,pi,pi) q[15];

// Unitary control gates
cx q[0], q[1];
cy q[2], q[3];
//cz q[0], q[1];
//ch q[1], q[2];
//csx q[0], q[1];

// Bell state
//h q[0];
//cx q[0], q[1];
//x q[0];

// Measure state

