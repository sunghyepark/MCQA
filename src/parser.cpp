//parser.cpp

#include "circuit.h"
#include "QASMtoken.hpp"
#include "QASMscanner.hpp"
#include "QASMparser.h"

#define PRINT_ParseOption 1

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::parsing(int argc, char** argv)
{
    if(argc != 5)
    {
        cout << "ERROR: invalid arguments" << endl;
        exit(1);
    }

    fileName_input  = argv[1];
    fileName_output = argv[2];
    
    string alpha = argv[3];
    string beta = argv[4];
    param_alpha = stod(alpha);
    param_beta = stod(beta);

#if PRINT_ParseOption
    cout << "*** Option ***" << endl;    
    cout << "Input file  : " << fileName_input << endl;    
    cout << "Output file  : " << fileName_output << endl;    
    cout << "parameter ( alpha: " << param_alpha << ", beta: " << param_beta << ")" << endl;
    cout << "" << endl;    
    cout << endl;    
#endif
    QASM_parse(fileName_input);
}

void Qcircuit::QMapper::QASM_parse(string filename)
{
    Circuit* circuit;
    // "*** Parsing Input File ***
    circuit = &Dgraph;
    
    //Start Parser
    QASMparser* parser = new QASMparser(filename);
    parser->Parse();
    nqubits = parser->getNqubits();
    circuit->nodeset = parser->getGatelists();
    //cout << circuit->nodeset << endl;
    delete parser;
}







