
#include "circuit.h"

#define PRINT_OUTPUT_info 1
#define FINALNODE_PRINT   0

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::FinalCircuit_info(Circuit& graph)
{
    //final circuit
#if PRINT_OUTPUT_info
    cout << "*** Final Circuit ***" << endl;
#if GATESET_IIC_JKU || GATESET_DAC2022
    cout << "# add_cnot_num: " << add_cnot_num << endl;
#elif GATESET_QUTECH
    cout << "# add_cz_num: " << add_cz_num << endl;
#endif
    cout << "# add_swap_num: " << add_swap_num << endl;
    cout << "# add_mov_num: "  << add_mov_num << endl;
    cout << "# add_bridge_num: " << add_bridge_num << endl;
    cout << "# latency: " << latency << endl;
#endif    

#if FINALNODE_PRINT
    cout << graph.nodeset << endl;
    cout << " ==> Print Final Circuit Successfully Finished" << endl;
#endif    
}    

void Qcircuit::QMapper::write_output(Circuit& graph)
{
    ofstream of(fileName_output);
    of << "OPENQASM 2.0;" << endl;
    of << "include \"qelib1.inc\";" << endl;
    of << "qreg q[20];" << endl; 
    of << "creg c[20];" << endl; 
    
    for (auto g : graph.nodeset)
    {
        if(g.control==-1)
            of << g.type <<  g.output_type << " " << "q[" << g.target << "];" << endl;
        else
            of << g.type << " " << "q[" << g.control << "],q[" << g.target << "];" << endl; 
    }
    of.close();
#if PRINT_OUTPUT_info
    cout << " ==> Output File (.qasm) Is Successfully Generated" << endl;
#endif

}

