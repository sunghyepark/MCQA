/*----------------------------------------------------------------------------------------------*/
/* Description:     Header file for QASM token                                                  */
/*                                                                                              */
/* Author:          Alwin Zulehner, Robert Wille                                                */
/*                                                                                              */
/* Revised by:      Sunghye Park, Postech                                                       */
/*                                                                                              */
/* Created:         02/10/2020 (for Quantum compiler)                                           */
/*----------------------------------------------------------------------------------------------*/
#ifndef TOKEN_H_
#define TOKEN_H_

#include <map>

#include "circuit.h"


class Token {
 public:

    enum class Kind {
        include, none, identifier, number, plus, semicolon, eof, lpar, rpar, lbrack, rbrack, 
        lbrace, rbrace, comma, minus, times, nninteger, real, qreg, creg,
#if GATESET_IIC_JKU
// >>>>>>>>>>>>>>>>>>>> IIC-JKU Original gateset >>>>>>>>>>>>>>>>>>>>
        ugate, cxgate,
// <<<<<<<<<<<<<<<<<<<< IIC-JKU Original gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_DAC2022
// >>>>>>>>>>>>>>>>>>>> DAC2022 target gateset >>>>>>>>>>>>>>>>>>>>
        xgate, ygate, zgate, hgate,
        sgate, sdggate, tgate, tdggate, 
        rxgate, rzgate, cxgate,
// <<<<<<<<<<<<<<<<<<<< DAC2022 target gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_QUTECH
// >>>>>>>>>>>>>>>>>>>> QuTech gateset >>>>>>>>>>>>>>>>>>>>
        rxgate, rygate, czgate,
// <<<<<<<<<<<<<<<<<<<< QuTech gateset <<<<<<<<<<<<<<<<<<<<
#endif // GATESET
        gate, pi, measure, openqasm, probabilities, sin, cos, tan, exp, 
        ln, sqrt, div, power, string, gt, barrier, opaque, _if, eq, reset, snapshot
    };

    Token(Kind kind, int line, int col) {
        this->kind = kind;
        this->line = line;
        this->col = col;
        this->val = 0;
        this->valReal = 0.0;
    }

    Token() : Token(Kind::none, 0, 0) {
    }

    static std::map<Kind, std::string> KindNames;
    Kind kind;
    int line;
    int col;
    int val;
    double valReal;
    std::string str;

 };

#endif /* TOKEN_H_ */
