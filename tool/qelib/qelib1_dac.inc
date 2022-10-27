// Quantum Experience (QE) Standard Header
// file: qelib1.inc

// controlled-NOT
gate cx c,t { CX c,t; }

// Pauli gate: bit-flip
gate x a { X a; }

// Pauli gate: bit and phase flip
gate y a { Y a; }

// Pauli gate: phase flip
gate z a { Z a; }

// Clifford gate: Hadamard
gate h a { H a; }

// Clifford gate: sqrt(Z) phase gate
gate s a { S a; }

// Clifford gate: conjugate of sqrt(Z)
gate sdg a { SDG a; }

// C3 gate: sqrt(Y) phase gate
gate t a { T a; }

// C3 gate: conjugate of sqrt(Y)
gate tdg a { TDG a; }
