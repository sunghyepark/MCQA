// Quantum Experience (QE) Standard Header
// file: qelib1.inc

// 3-parameter 2-pulse single qubit gate
gate u3(theta,phi,lambda) q { U(theta,phi,lambda) q; }

// 2-parameter 1-pulse single qubit gate
gate u2(phi,lambda) q { U(pi/2,phi,lambda) q; }

// 1-parameter 0-pulse single qubit gate
gate u1(lambda) q { U(0,0,lambda) q; }

// controlled-NOT
gate cx c,t { CX c,t; }

// idle gate (identity)
gate id a { U(0,0,0) a; }

// Pauli gate: bit-flip
gate x a { u3(pi,0,pi) a; }

// Pauli gate: bit and phase flip
gate y a { u3(pi,pi/2,pi/2) a; }

// Pauli gate: phase flip
gate z a { u1(pi) a; }

// Clifford gate: Hadamard
gate h a { u2(0,pi) a; }

// Clifford gate: sqrt(Z) phase gate
gate s a { u1(pi/2) a; }

// Clifford gate: conjugate of sqrt(Z)
gate sdg a { u1(-pi/2) a; }

// C3 gate: sqrt(Y) phase gate
gate t a { u1(pi/4) a; }

// C3 gate: conjugate of sqrt(Y)
gate tdg a { u1(-pi/4) a; }
