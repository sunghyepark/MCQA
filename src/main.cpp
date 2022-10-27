// main.cpp

/* --------------------------------------- */
/* This script is written by Sunghye Park  */
/* 2022.03.21                              */
/*                                         */
/* shpark96@postech.ac.kr                  */
/* --------------------------------------- */

#include "circuit.h"
#include "QASMtoken.hpp"
#include "QASMscanner.hpp"
#include "QASMparser.h"
#include "mymeasure.h"

    //execute each procedure
#define MCQC_initial 1
#define MCQC_main    1
#define MCQC_out     1
#define MCQC_pre     1  // dhkim, 2022.05.03
#define MCQC_post    1  // dhkim, 2022.05.03

    //coupling_graph
#define ARCHI ARCHITECTURE::Surface_17
    //initial mapping
#define INITIAL INITIALTYPE::GRAPH_MATCHING_PROPOSED
#define PRE_PROCESSING 1
    //main mapping
#define BRIDGE_MODE    1

using namespace std;
using namespace Qcircuit;

Qcircuit::QMapper mapper;
static CMeasure measure;

int main(int argc, char** argv)
{
    cout << endl;
    cout << "======================================================================================" << endl;
    cout << "             MCQA (multi-constraint quantum allocation) by Sunghye Park               " << endl;
    cout << "                                                           Dohun   Kim                " << endl;
    cout << "======================================================================================" << endl;
    cout << endl;
    measure.start_clock();
    
    //parsing
    cout << "\n\t================== 1. Circuit parsing ========================" << endl;
    mapper.parsing(argc, argv);
    measure.stop_clock("QASM parsing > graph");

#if MCQC_pre
    // preprocessing: single-qubit gate optimization
    // NOTE: added by dhkim, 2022.04.28
    mapper.circuit_optimization(mapper.Dgraph, true);
    measure.stop_clock("Preprocessing");
#endif

    //mapping
    cout << "\n\t================== 2. Initial mapping ========================" << endl;
#if MCQC_initial
    mapper.initial_mapping(ARCHI, PRE_PROCESSING, BRIDGE_MODE, INITIAL);
    measure.stop_clock("Initial Mapping");
#endif

    cout << "\n\t================== 3. Main mapping ============================" << endl;
#if MCQC_main
    mapper.Circuit_Mapping(mapper.Dgraph, ARCHI, BRIDGE_MODE);
    measure.stop_clock("Main Mapping");
#endif

   cout << "\n\t================== 3-2. Post optimization ========================" << endl;
#if MCQC_post
    // postprocessing: single-qubit gate optimization
    // NOTE: added by dhkim, 2022.05.03
    mapper.circuit_optimization(mapper.FinalCircuit, false);
    mapper.update_latency(mapper.FinalCircuit, ARCHI);    
    measure.stop_clock("Post optimization");
#endif

    cout << "\n\t================== 4. Output write ============================" << endl;
#if MCQC_out
    mapper.FinalCircuit_info(mapper.FinalCircuit);
    mapper.write_output(mapper.FinalCircuit);
    measure.stop_clock("Print Output");
#endif

    //measure
    measure.stop_clock("END Program");
    cout << "\n*** Time & Memory Measure ***" << endl;
    measure.print_clock();
    measure.printMemoryUsage();

    return 0;
}
