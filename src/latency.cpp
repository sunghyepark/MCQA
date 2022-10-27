//layout.cpp
//added by dhkim, 2022.05.03

#include "circuit.h"

namespace Qcircuit {

void QMapper::update_latency(Circuit& dgraph, ARCHITECTURE archi) {

    struct result_values {
#if GATESET_IIC_JKU || GATESET_DAC2022
        int n_cnot;
#elif GATESET_QUTECH
        int n_cz;
#endif
        int n_swap;
        int n_mov;
        int n_bridge;
    };

    result_values val{
#if GATESET_IIC_JKU || GATESET_DAC2022
        add_cnot_num,
#elif GATESET_QUTECH
        add_cz_num,
#endif
        add_swap_num,
        add_mov_num,
        add_bridge_num
    };

    Dgraph = dgraph;
    layout = initial_layout;

    Circuit_Mapping(Dgraph, archi, 1);

#if GATESET_IIC_JKU || GATESET_DAC2022
        add_cnot_num = val.n_cnot;
#elif GATESET_QUTECH
        add_cz_num = val.n_cz;
#endif
        add_swap_num = val.n_swap;
        add_mov_num = val.n_mov;
        add_bridge_num = val.n_bridge;
}

} // namespace Qcircuit
