// QuTech Primitive Gateset Header
// file: qelib1.inc

// --- QuTech Hardware primitives ---

gate rx45 a { RX(pi/4) a; }
gate rx90 a { RX(pi/2) a; }
gate rx180 a { RX(pi) a; }

gate rxm45 a { RX(-pi/4) a; }
gate rxm90 a { RX(-pi/2) a; }
gate rxm180 a { RX(-pi) a; }

gate ry45 a { RY(pi/4) a; }
gate ry90 a { RY(pi/2) a; }
gate ry180 a { RY(pi) a; }

gate rym45 a { RY(-pi/4) a; }
gate rym90 a { RY(-pi/2) a; }
gate rym180 a { RY(-pi) a; }

gate cz c,t { CZ c,t; }

// --- QuTech Standard Gates ---

gate x a { rx180 a; }
gate y a { ry180 a; }
gate z a { rx180 a; ry180 a; }
gate h a { ry90 a; }
gate cx a,b { rym90 b; cz a,b; ry90 b; }

gate s a { ry90 a; rx90 a; rym90 a; }
gate sdg a { ry90 a; rxm90 a; rym90 a; }

gate t a { ry90 a; rx45 a; rym90 a; }
gate tdg a { ry90 a; rxm45 a; rym90 a; }


