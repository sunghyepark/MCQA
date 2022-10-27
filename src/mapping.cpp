//mapping.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

bool compare_pairSecond( pair<int, int> p1, pair<int, int> p2 )
{
    return p1.second < p2.second;
}

void Qcircuit::QMapper::Circuit_Mapping(Circuit& dgraph, ARCHITECTURE archi, bool BRIDGE_MODE)
{
    //////////// (0) make Circuit graph ////////////
    make_Glist(dgraph);

    //// (0-1) make frequency relationship
    vector<pair<vector<int>, vector<int>>>  freq_relations; //[m]: physical qubits -> first: lower freq.groups/ second: higher freq. groups
    vector<vector<int>> same_freq_qubits;               //[m]: same_freq_qubits
    if(archi == ARCHITECTURE::Surface_17)
    {
        generate_freq_relations(freq_relations);
        generate_freq_same(same_freq_qubits);
    }

    //////////// (1) Circuit mapping ///////////////
    //cost weight table -----------------
    cost_weight_table.clear();
    double cost = 1.0;

    //initiallize ------------------------
        //gate list
    list<int> g_front;
    list<int> g_pair;
    list<int> g_sat;
    list<int> g_vio;
    list<int> g_bridge;
    list<int> g_scheduled;
    vector<pair<pair<int, int>, pair<int,int>> > g_temp; //<gateid, qubits>, <two or single, length of Glist>>

        //frozen
    vector<int> frozen(positions, 0);                  //0: not frozen --> update gate is possible, 1: frozen --> update gate is impossible
    vector<int> frozen_duration(positions, 0);    //1: single-qubit gate, 2: two-qubit gate
    vector<int> frozen_frequency(positions, -1);  //0: original, 1: interaction, 2:parking, 3: single-qubit
    vector<int> frozen_constraints(positions, 0);       //0: not frozen, 1: frozen
        
        //iteration
    int loop_end = 0;

        //output circuit
    FinalCircuit.nodeset.clear();
        
        //resolve constraint
    node_id = dgraph.nodeset.size();    //for add_SWAP
#if GATESET_IIC_JKU || GATESET_DAC2022
    add_cnot_num   = 0;
#elif GATESET_QUTECH
    add_cz_num     = 0;
#endif
    add_swap_num   = 0;
    add_mov_num    = 0;
    add_bridge_num = 0;
    latency = 0;
        
    vector< pair< pair<int,int>, pair<int, double> > > MCPE_swap_flag;   // SWAP candidate, act_gate_id, cost
    vector< pair< pair<int,int>, pair<int, double> > > MCPE_mov_flag;    // SWAP candidate, act_gate_id, cost
    int history_size = 2;

        //for layout update
    vector< pair<int, int> > scheduler_constraints; //(layout_swap_n_mov, number_swap_n_mov) //0:BRIDGE/ 1:SWAP/ 2:MOV //  
    vector< pair<pair<int,int>, int> > qubit_constraints;     //(q_c, q_t),q_m
    vector<bool> updated_Glist(nqubits, false);
    vector< pair<pair<int,int>, int> > gateinfo_constraints; //gateid for q_c, q_t, q_m
    vector< pair<pair<int,int>, int> > gateexecute_constraints; //execute final gateid

    //do while front_list[] is not empty -----------------------------
    do{
        //cout << "\n----------------- latency: " << latency+1 << " -----------------" << endl;
        for(int i=0; i<nqubits; i++)
            updated_Glist[i] = false;

        update_frozen_duration(frozen_duration);
        update_frozen_frequency(frozen_frequency);
       
        //(1)  prepare gate lists ---------------------
        // (1-1) Update front & pair gates
        update_front_n_pair_gates(g_front, g_pair, g_sat, g_vio, frozen);
        
        // (1-2) Update sat & vio gates
        bool complete_g_pair = true;
        complete_g_pair = check_connectivity(g_front, g_pair, g_sat, g_vio, frozen, dgraph);

        //(2) constraint not exist -> push ------------
        if(!g_sat.empty())
        {
            //print_gates(g_sat, 2);
        }
        else if(!g_vio.empty())
        {
            // (3-0-1) check g_vio
            latency--;

            //(3-0-2) scheduling
            generate_scheduled_duration(g_vio, g_temp, frozen_duration, dgraph);
            sort_priority_gates(g_temp, g_scheduled, frozen_duration);
            //(3-1) make candidates ----------------------
                //calculate qubit priority
            vector< pair<int, int> > qubit_length_pair;
            vector<int> qubit_length_priority;
            for(int n=0; n<nqubits; n++)
            {
                int _Dlength = Glist[n].size();
                if(_Dlength != 0)
                    qubit_length_pair.push_back( make_pair(n, _Dlength) );
            }
            if(!qubit_length_pair.empty())
            {
                sort(qubit_length_pair.begin(), qubit_length_pair.end(), compare_pairSecond);
                for(auto& n : qubit_length_pair)
                    qubit_length_priority.push_back(n.first);
            }

                //3-1 SWAP //////////////
            vector< pair< pair<int, int>, int> > swap_candi_list; //SWAP candidates, g_vio id
            generate_swap_candi_list(g_vio, swap_candi_list, frozen_constraints, dgraph);
            
            //MCPE_cost
            vector< pair< pair<int, int>, pair<int, double> > > MCPE_swap_test; //SWAP candidate, g_vio_id, cost
            for(auto kv : swap_candi_list)
            {
                pair<int, int> SWAP_pair = kv.first;
                double cost = cal_MCPE(SWAP_pair, dgraph, true);
                MCPE_swap_test.push_back(make_pair(SWAP_pair, make_pair(kv.second, cost)));
                //cout << "MCPE_swap_test: " << MCPE_swap_test << endl;
            }
                
                //3-2 MOV ///////////////
            vector< pair< pair<int, int>, int> > mov_candi_list; //MOV candidates, g_vio id
            generate_mov_candi_list(g_vio, mov_candi_list, frozen_constraints, dgraph);
                
            //MCPE_cost
            vector< pair< pair<int, int>, pair<int, double> > > MCPE_mov_test; //MOV candidate, g_vio_id, cost
            for(auto kv : mov_candi_list)
            {
                pair<int, int> MOV_pair = kv.first;
                double cost = cal_MCPE(MOV_pair, dgraph, false);
                MCPE_mov_test.push_back(make_pair(MOV_pair, make_pair(kv.second, cost)));
                //cout << "MCPE_mov_test: " << MCPE_mov_test << endl;
            }
                
                    //3-3 BRIDGE ////////////
            if(BRIDGE_MODE)
            {
                update_bridge_gates(g_bridge, g_vio, frozen_constraints, dgraph);
                vector<pair<pair<int, int>, pair<int,int>> > temp_gates; //<gateid, qubits>, <two or single, length of Glist>>
                generate_scheduled_duration(g_bridge, temp_gates, frozen_duration, dgraph);
                sort_priority_gates(temp_gates, g_bridge, frozen_duration);
            }
            
            //(3-2) constraint resolved ----------------------
            pair<int, int> SWAP;
            int swap_gateid;
            double max_swap_cost = 0.0;
            if(!MCPE_swap_test.empty())
            {
                find_max_cost(SWAP, swap_gateid, max_swap_cost, MCPE_swap_test, g_vio, MCPE_swap_flag);
                //print_candi_pair(SWAP, swap_gateid, max_swap_cost, 0);
            }

            pair<int, int> MOV;
            int mov_gateid;
            double max_mov_cost = 0.0;
            if(!MCPE_mov_test.empty())
            {
                find_max_cost(MOV, mov_gateid, max_mov_cost, MCPE_mov_test, g_vio, MCPE_mov_flag);
                //print_candi_pair(MOV, mov_gateid, max_mov_cost, 0);
            }
            
            //max cost?? (between swap and mov)
            bool is_swap = false;
            bool is_mov  = false;
            int gateid;
            double max_cost = 0.0;
            if(max_swap_cost > max_mov_cost && max_swap_cost > 0.0)
            {
                //cout << "max_swap_cost: " << max_swap_cost << endl;
                is_swap  = true;
                gateid   = swap_gateid;
                max_cost = max_swap_cost;
            }
            else if(max_mov_cost > 0.0)
            {
                //cout << "max_mov_cost: " << max_mov_cost << endl;
                is_mov   = true;
                gateid   = mov_gateid;
                max_cost = max_mov_cost;
            }
            else
                max_cost = 0.0;
            
            //resolve constraint //////////////////////////////////////
            if(max_cost == 0.0)
                gateid = g_bridge.front();
            std::list<int>::iterator it = std::find(g_bridge.begin(), g_bridge.end(), gateid);
            if(it != g_bridge.end() && max_cost <= param_beta * 10) 
            {
                //3-3 BRIDGE ////////////
                const int q_B_control = dgraph.nodeset[gateid].control;
                const int q_B_target  = dgraph.nodeset[gateid].target;
                const int p_control = layout.log2phy(q_B_control);
                const int p_target  = layout.log2phy(q_B_target );
                //cout << "\n -> add_bridge" << endl;
                //print_addGates(gateid, q_B_control, q_B_target, p_control, p_target);
                
                //remove @Glist
                Glist[q_B_control].pop_front();
                Glist[q_B_control].pop_front();
                Glist[q_B_target ].pop_front();
                Glist[q_B_target ].pop_front();
                
                frozen[q_B_control] = 0;
                frozen[q_B_target ] = 0;

                //add_bridge_Glist();
                add_bridge(p_control, p_target, frozen_constraints, dgraph);
#if GATESET_IIC_JKU || GATESET_DAC2022
                gateinfo_constraints.push_back(make_pair( make_pair( node_id-1, node_id-2), node_id-1) );
#elif GATESET_QUTECH
                gateinfo_constraints.push_back(make_pair( make_pair( node_id-2, node_id-4), node_id-1) );
#endif
                gateexecute_constraints.push_back(make_pair( make_pair( 0, 0 ), 0 ) );
                int bridge_nodeid = node_id;

                const int q_B_middle = dgraph.nodeset[node_id-1].target;
                //cout << "q_B_middle" << q_B_middle << endl;
                int bridge_related_gateid = Glist[q_B_middle].front().first;

                frozen[q_B_middle ] = 0;
                
                //remove original gateid
                g_vio.remove(gateid);
                g_pair.remove(gateid);
                
                //remove bridge_related gateid
                if(Glist[q_B_middle].size() != 0)
                    g_front.remove(bridge_related_gateid);
                
                //Glist update (descending order!!)
                add_gate_to_Glist(0, bridge_nodeid, q_B_control, q_B_target, q_B_middle);
                scheduler_constraints.push_back( make_pair(0, 1) ); //layout_swap_n_mov & num_swap_n_mov
                qubit_constraints.push_back( make_pair( make_pair(q_B_control, q_B_target), q_B_middle));
                frozen_constraints[p_control] = 1;
                frozen_constraints[p_target ] = 1;
                const int p_B_middle = layout.log2phy(q_B_middle);
                frozen_constraints[ p_B_middle ] = 1;
                //cout << "q_control: " << q_B_control << ", q_target: " << q_B_target << ", q_middle: " << q_B_middle << endl;
                //cout << "p_control: " << p_control   << ", p_target: " << p_target   << ", p_middle: " << p_B_middle << endl;
            }
            else
            {
                if(is_swap)
                {
                    //3-1 SWAP //////////////
                    int q_control = dgraph.nodeset[gateid].control;
                    int q_target  = dgraph.nodeset[gateid].target;
                    int p_control = layout.log2phy(q_control);
                    int p_target  = layout.log2phy(q_target );
                    //cout << "\n -> add_swap: " << endl;
                    //print_addGates(gateid, q_control, q_target, p_control, p_target);
                    
                    //remove @Glist
                    frozen[q_control] = 0;
                    frozen[q_target ] = 0;
                    
                    
                    //add swap
                    if(MCPE_swap_flag.size() > history_size)
                        MCPE_swap_flag.erase(MCPE_swap_flag.begin());
                    MCPE_swap_flag.push_back(MCPE_swap_test[0]);

                    //cout << "** we SWAP " << SWAP << endl;
                    int q_SWAP_first  = layout.phy2log(SWAP.first );
                    int q_SWAP_second = layout.phy2log(SWAP.second);

                    add_swap(q_SWAP_first, q_SWAP_second, dgraph);
#if GATESET_IIC_JKU || GATESET_DAC2022
                    gateinfo_constraints.push_back(make_pair( make_pair( node_id-1, node_id-1), -1) );
#elif GATESET_QUTECH
                    gateinfo_constraints.push_back(make_pair( make_pair( node_id-2, node_id-1), -1) );
#endif
                    gateexecute_constraints.push_back(make_pair( make_pair( 0, 0 ), 0 ) );

                    int swap_nodeid = node_id;

                    //remove original gateid
                    g_vio.remove(gateid);
                    g_pair.remove(gateid);
                    
                    //remove swap_related gateid
                    int q_middle;
                    int swap_related_gateid;
                    if(q_control == q_SWAP_first)
                        q_middle = q_SWAP_second;
                    else if(q_target == q_SWAP_first)
                        q_middle = q_SWAP_second;
                    else if(q_control == q_SWAP_second)
                        q_middle = q_SWAP_first;
                    else if(q_target == q_SWAP_second)
                        q_middle = q_SWAP_first;
                    swap_related_gateid = Glist[q_middle].front().first;
                    
                    if(Glist[q_middle].size() != 0)
                    {
                        g_front.remove(swap_related_gateid);
                        
                        std::list<int>::iterator it = std::find(g_pair.begin(), g_pair.end(), swap_related_gateid);
                        if(it != g_pair.end()) 
                        {
                            g_pair.remove(swap_related_gateid);
                            int q_control_swap = Dgraph.nodeset[swap_related_gateid].control;
                            int q_target_swap  = Dgraph.nodeset[swap_related_gateid].target ;
                            frozen[q_control_swap] = 0;
                            frozen[q_target_swap ] = 0;

                            std::list<int>::iterator it2 = std::find(g_vio.begin(), g_vio.end(), swap_related_gateid);
                            if(it2 != g_vio.end()) 
                                g_vio.remove(swap_related_gateid);
                        }
                    }
                    frozen[q_middle ] = 0;
                    
                    //Glist update (descending order!!)
                    add_gate_to_Glist(1, swap_nodeid, q_SWAP_first, q_SWAP_second, 0);
                    
                    //layout swap
                    scheduler_constraints.push_back( make_pair(1, 1) ); //layout_swap_n_mov & num_swap_n_mov
                    qubit_constraints.push_back( make_pair( make_pair(q_SWAP_first, q_SWAP_second), -1));
                    
                    frozen_constraints[SWAP.first]  = 1;
                    frozen_constraints[SWAP.second] = 1;
                }
                if(is_mov)
                {
                    //3-2 MOV ///////////////
                    const int q_control = dgraph.nodeset[gateid].control;
                    const int q_target  = dgraph.nodeset[gateid].target;
                    const int p_control = layout.log2phy(q_control);
                    const int p_target  = layout.log2phy(q_target );
                    //cout << "\n -> add_mov: " << endl;
                    //print_addGates(gateid, q_control, q_target, p_control, p_target);
                    
                    //remove @Glist
                    frozen[q_control] = 0;
                    frozen[q_target ] = 0;
                    
                    //add mov
                    if(MCPE_mov_flag.size() > history_size)
                        MCPE_mov_flag.erase(MCPE_mov_flag.begin());
                    MCPE_mov_flag.push_back(MCPE_mov_test[0]);

                    //cout << "** we MOV " << MOV << endl;
                    const int q_MOV_first  = layout.phy2log(MOV.first );
                    const int q_MOV_second = layout.phy2log(MOV.second);
                    add_mov(q_MOV_first, q_MOV_second, dgraph);
#if GATESET_IIC_JKU || GATESET_DAC2022
                    gateinfo_constraints.push_back(make_pair( make_pair( node_id-1, node_id-1), -1) );
#elif GATESET_QUTECH
                    gateinfo_constraints.push_back(make_pair( make_pair( node_id-1, node_id-2), -1) );
#endif
                    gateexecute_constraints.push_back(make_pair( make_pair( 0, 0 ), 0 ) );

                    int mov_nodeid = node_id;
                    
                    //remove original gateid
                    g_vio.remove(gateid);
                    g_pair.remove(gateid);

                    //Glist update (descending order!!)
                    add_gate_to_Glist(2, mov_nodeid, q_MOV_first, q_MOV_second, 0);
                    
                    //layout mov
                    scheduler_constraints.push_back( make_pair(2, 1) ); //layout_swap_n_mov & num_swap_n_mov
                    qubit_constraints.push_back( make_pair( make_pair(q_MOV_first, q_MOV_second), -1));
                    
                    frozen_constraints[MOV.first]  = 1;
                    frozen_constraints[MOV.second] = 1;
                }
            }
        }
            
       //(4-1) scheduling??? -----
        generate_scheduled_duration(g_sat, g_temp, frozen_duration, dgraph);
        sort_priority_gates(g_temp, g_scheduled, frozen_duration);
        generate_scheduled_frequency(g_scheduled, frozen_frequency, dgraph);
        //(4-2) push to final circuit
        directly_execute_gate(g_scheduled, g_sat, frozen, frozen_duration, updated_Glist, dgraph);

        // (4-3) After executing gate
        //print_gates(g_scheduled, 5);
        update_Glist(frozen_duration, updated_Glist, dgraph);

        //(5) scheduling swap & mov
        vector<int> idx_schedule;
        vector<int> pop_idx;
        idx_schedule.clear();
        pop_idx.clear();
        for(int i=0; i<qubit_constraints.size(); i++)
        {
            const int q_c = qubit_constraints[i].first.first ;
            const int q_t = qubit_constraints[i].first.second;
#if GATESET_IIC_JKU || GATESET_DAC2022
            if(scheduler_constraints[i].first != 0)
            {
                //cout << "swap or mov" << endl;
                if(!gateexecute_constraints[i].first.first)
                {
                    int gate_c = gateinfo_constraints[i].first.first ;
                    list<int>::iterator it_c = find(g_scheduled.begin(), g_scheduled.end(), gate_c);
                    if(it_c != g_scheduled.end())
                        gateexecute_constraints[i].first.first  = 1;
                }
                else if(gateexecute_constraints[i].first.first == 1 && updated_Glist[q_c] == true)
                {
                    gateexecute_constraints[i].first.first = 2;
                    idx_schedule.push_back(i);
                }
                if(!gateexecute_constraints[i].first.second)
                {
                    int gate_t = gateinfo_constraints[i].first.second;
                    list<int>::iterator it_t = find(g_scheduled.begin(), g_scheduled.end(), gate_t);
                    if(it_t != g_scheduled.end())
                        gateexecute_constraints[i].first.second = 1;
                }
                else if(gateexecute_constraints[i].first.second == 1 && updated_Glist[q_t] == true)
                {
                    gateexecute_constraints[i].first.second = 2;
                    idx_schedule.push_back(i);
                }
            }
            else
            {
                //cout << "bridge" << endl;
                const int q_m = qubit_constraints[i].second;
                if(!gateexecute_constraints[i].first.first)
                {
                    int gate_c = gateinfo_constraints[i].first.first ;
                    list<int>::iterator it_c = find(g_scheduled.begin(), g_scheduled.end(), gate_c);
                    if(it_c != g_scheduled.end())
                        gateexecute_constraints[i].first.first  = 1;
                }
                else if(gateexecute_constraints[i].first.first == 1 && updated_Glist[q_c] == true)
                {
                    gateexecute_constraints[i].first.first = 2;
                    idx_schedule.push_back(i);
                }
                if(!gateexecute_constraints[i].first.second)
                {
                    int gate_t = gateinfo_constraints[i].first.second;
                    list<int>::iterator it_t = find(g_scheduled.begin(), g_scheduled.end(), gate_t);
                    if(it_t != g_scheduled.end())
                        gateexecute_constraints[i].first.second = 1;
                }
                else if(gateexecute_constraints[i].first.second == 1 && updated_Glist[q_t] == true)
                {
                    gateexecute_constraints[i].first.second = 2;
                    idx_schedule.push_back(i);
                }
                if(!gateexecute_constraints[i].second)
                {
                    int gate_m = gateinfo_constraints[i].second;
                    list<int>::iterator it_m = find(g_scheduled.begin(), g_scheduled.end(), gate_m);
                    if(it_m != g_scheduled.end())
                        gateexecute_constraints[i].second = 1;
                }
                else if(gateexecute_constraints[i].second == 1 && updated_Glist[q_m] == true)
                {
                    gateexecute_constraints[i].second = 2;
                    idx_schedule.push_back(i);
                }
            }
#elif GATESET_QUTECH
            if(scheduler_constraints[i].first == 1)
            {
                //cout << "swap" << endl;
                if(!gateexecute_constraints[i].first.first)
                {
                    int gate_c = gateinfo_constraints[i].first.first ;
                    list<int>::iterator it_c = find(g_scheduled.begin(), g_scheduled.end(), gate_c);
                    if(it_c != g_scheduled.end())
                        gateexecute_constraints[i].first.first  = 1;
                }
                else if(gateexecute_constraints[i].first.first == 1 && updated_Glist[q_c] == true)
                {
                    gateexecute_constraints[i].first.first = 2;
                    idx_schedule.push_back(i);
                }
                if(!gateexecute_constraints[i].first.second)
                {
                    int gate_t = gateinfo_constraints[i].first.second;
                    list<int>::iterator it_t = find(g_scheduled.begin(), g_scheduled.end(), gate_t);
                    if(it_t != g_scheduled.end())
                    {
                        gateexecute_constraints[i].first.second = 1;
                        idx_schedule.push_back(i);
                    }
                }
            }
            if(scheduler_constraints[i].first == 2)
            {
                //cout << "mov" << endl;
                if(!gateexecute_constraints[i].first.first)
                {
                    int gate_c = gateinfo_constraints[i].first.first ;
                    list<int>::iterator it_c = find(g_scheduled.begin(), g_scheduled.end(), gate_c);
                    if(it_c != g_scheduled.end())
                    {
                        gateexecute_constraints[i].first.first  = 1;
                        idx_schedule.push_back(i);
                    }
                }
                if(!gateexecute_constraints[i].first.second)
                {
                    int gate_t = gateinfo_constraints[i].first.second;
                    list<int>::iterator it_t = find(g_scheduled.begin(), g_scheduled.end(), gate_t);
                    if(it_t != g_scheduled.end())
                        gateexecute_constraints[i].first.second = 1;
                }
                else if(gateexecute_constraints[i].first.second == 1 && updated_Glist[q_t] == true)
                {
                    gateexecute_constraints[i].first.second = 2;
                    idx_schedule.push_back(i);
                }
            }
            else
            {
                //cout << "bridge" << endl;
                const int q_m = qubit_constraints[i].second;
                if(!gateexecute_constraints[i].first.first)
                {
                    int gate_c = gateinfo_constraints[i].first.first ;
                    list<int>::iterator it_c = find(g_scheduled.begin(), g_scheduled.end(), gate_c);
                    if(it_c != g_scheduled.end())
                        gateexecute_constraints[i].first.first  = 1;
                }
                else if(gateexecute_constraints[i].first.first == 1 && updated_Glist[q_c] == true)
                {
                    gateexecute_constraints[i].first.first = 2;
                    idx_schedule.push_back(i);
                }
                if(!gateexecute_constraints[i].first.second)
                {
                    int gate_t = gateinfo_constraints[i].first.second;
                    list<int>::iterator it_t = find(g_scheduled.begin(), g_scheduled.end(), gate_t);
                    if(it_t != g_scheduled.end())
                    {
                        gateexecute_constraints[i].first.second = 1;
                        idx_schedule.push_back(i);
                    }
                }
                if(!gateexecute_constraints[i].second)
                {
                    int gate_m = gateinfo_constraints[i].second;
                    list<int>::iterator it_m = find(g_scheduled.begin(), g_scheduled.end(), gate_m);
                    if(it_m != g_scheduled.end())
                    {
                        gateexecute_constraints[i].second = 1;
                        idx_schedule.push_back(i);
                    }
                }
            }
#endif
        }
        if(!idx_schedule.empty())
        {
            for(auto& i : idx_schedule)
            {
                //layout_swap_n_mov == 0 (BRIDGE)
                if(scheduler_constraints[i].first == 0)
                {
                    int q_control = qubit_constraints[i].first.first;
                    int q_target  = qubit_constraints[i].first.second;
                    int q_middle  = qubit_constraints[i].second;
                    //cout << "* bridge (" << q_control << ", " << q_target << ") -> " << scheduler_constraints[i].second << endl; 
                    if(scheduler_constraints[i].second == 3)
                    {
                        int p_B_control = layout.log2phy(q_control);
                        int p_B_target  = layout.log2phy(q_target );
                        int p_B_middle  = layout.log2phy(q_middle ); 
                        //cout << "q_control: " << q_control   << ", q_target: " << q_target   << ", q_middle: " << q_middle   << endl;
                        //cout << "p_control: " << p_B_control << ", p_target: " << p_B_target << ", p_middle: " << p_B_middle << endl;
                        frozen_constraints[p_B_control] = 0;
                        frozen_constraints[p_B_target ] = 0;
                        frozen_constraints[p_B_middle ] = 0;
                        //cout << "* bridge is done" << endl;
                        pop_idx.push_back(i);
                    }
                    scheduler_constraints[i].second += 1;
                }
                //layout_swap_n_mov == 1 (SWAP)
                if(scheduler_constraints[i].first == 1) 
                {
                    int q_first  = qubit_constraints[i].first.first ;
                    int q_second = qubit_constraints[i].first.second;
                    //cout << "* swap (" << q_first << ", " << q_second << ") -> " << scheduler_constraints[i].second << endl; 
                    if(scheduler_constraints[i].second == 2)
                    {
                        layout.swap_logical(q_first, q_second);
                        int p_SWAP_first  = layout.log2phy(q_first );
                        int p_SWAP_second = layout.log2phy(q_second);
                        frozen_constraints[p_SWAP_first ] = 0;
                        frozen_constraints[p_SWAP_second] = 0;
                        //cout << "* layout swap" << endl;
                        //layout.print_log2phy();
                        pop_idx.push_back(i);
                    }
                    scheduler_constraints[i].second += 1; 
                }
                //layout_swap_n_mov == 2 (MOV)
                if(scheduler_constraints[i].first == 2)
                {
                    int q_first  = qubit_constraints[i].first.first ;
                    int q_second = qubit_constraints[i].first.second;
                    //cout << "* mov (" << q_first << ", " << q_second << ") -> " << scheduler_constraints[i].second << endl; 
                    if(scheduler_constraints[i].second == 2)
                    {
                        layout.swap_logical(q_first, q_second);
                        int p_MOV_first  = layout.log2phy(q_first );
                        int p_MOV_second = layout.log2phy(q_second);
                        frozen_constraints[p_MOV_first ] = 0;
                        frozen_constraints[p_MOV_second] = 0;
                        //cout << "* layout mov" << endl;
                        //layout.print_log2phy();
                        pop_idx.push_back(i);
                    }
                    scheduler_constraints[i].second += 1; 
                }
            }
            //
            if(!pop_idx.empty())
            {
                reverse(pop_idx.begin(), pop_idx.end());
                for(auto& i : pop_idx)
                {
                    qubit_constraints.erase( qubit_constraints.begin() + i );
                    scheduler_constraints.erase( scheduler_constraints.begin() + i );
                    gateinfo_constraints.erase( gateinfo_constraints.begin() + i );
                    gateexecute_constraints.erase( gateexecute_constraints.begin() + i );
                }
            }
        }
        
        //(6) for next iteration
        loop_end = 0;
        for(int q=0; q<nqubits; q++)
            if(Glist[q].empty()) loop_end++;
        latency++;
    }while(true && (loop_end!=nqubits) );
}

void Qcircuit::QMapper::directly_execute_gate(list<int>& g_scheduled, list<int>& g_sat, vector<int>& frozen, vector<int>& frozen_duration, vector<bool>& updated_Glist, Circuit& dgraph)
{
    
    // NOTE: added by dhkim, 2022.05.03
    // finds initial logical qubit: updated_layout -> initial_layout
    auto to_initial_qubit = [&](int q)
    {
        return initial_layout.phy2log(layout.log2phy(q));
    };

    vector<int> executed_gates;
    for(auto& gateid : g_scheduled)
    {
        const int q_control = dgraph.nodeset[gateid].control;
        const int q_target  = dgraph.nodeset[gateid].target;
        const int p_control = layout.log2phy(q_control);
        const int p_target  = layout.log2phy(q_target );
        const int f_control = coupling_graph.nodeset[p_control].getweight();
        const int f_target  = coupling_graph.nodeset[p_target ].getweight();
        if(q_control != -1) //two-qubit gate
        {
            //push to final circuit
            Gate gate = dgraph.nodeset[gateid];
            gate.control = to_initial_qubit(gate.control);
            gate.target = to_initial_qubit(gate.target);
            FinalCircuit.nodeset.push_back(gate);
            executed_gates.push_back(gateid);
            
            //update frozen gate_duration & frozen gate
            frozen_duration[q_control] = 2;
            frozen_duration[q_target]  = 2;
            frozen[q_control] = 0;
            frozen[q_target]  = 0;

            ////remove @g_sat
            g_sat.remove(gateid);

            //remove @Glist
            Glist[q_control].pop_front();
            Glist[q_target ].pop_front();

            //update
            updated_Glist[q_control]  = true;
            updated_Glist[q_target ]  = true;
        }
        else                //single-qubit gate
        {
            //push to final circuit
            Gate gate = dgraph.nodeset[gateid];
            gate.target = to_initial_qubit(gate.target);
            FinalCircuit.nodeset.push_back(gate);
            executed_gates.push_back(gateid);
            
            //update frozen gate_duration & frozen gate
            frozen_duration[q_target]  = 1;
            frozen[q_target]  = 0;

            ////remove @g_sat
            g_sat.remove(gateid);

            //remove @Glist
            Glist[q_target ].pop_front();

            //update
            updated_Glist[q_target ] = true;
        }
    }
}

void Qcircuit::QMapper::update_Glist(vector<int>& frozen_duration, vector<bool>& updated_Glist, Circuit& dgraph)
{
    for(int q=0; q<positions; q++)
    {   
        if (Glist[q].empty())
            continue;
        int gateid   = Glist[q].front().first;
        int frozen_q = frozen_duration[q];
        if(frozen_q == 1)
        {
            if(gateid == -1)  //temp of two-qubit gate
            {
                Glist[q].pop_front();
                updated_Glist[q]  = true;
            }
        }
    }
}

//frontgates, pairgates, frozen
void Qcircuit::QMapper::update_front_n_pair_gates(list<int>& g_front, list<int>& g_pair, list<int>& g_sat, list<int>& g_vio, vector<int>& frozen)
{
    for(int q=0; q<positions; q++)
    {
        list<pair<int, int>>& Glist_line = Glist[q];
        if(!frozen[q])
        {
            if(Glist_line.empty()) continue;
            int gateid = Glist_line.front().first;
            // 1. gateid = -1 -> pop
            if(gateid != -1)
            {
                // 2. update front & pair gates
                if( find(g_front.begin(), g_front.end(), gateid) != g_front.end() )
                {
                    g_front.remove(gateid);
                    g_pair.push_back(gateid);
                }
                else
                    g_front.push_back(gateid);
                frozen[q] = 1;
            }
        }
    }
}

bool Qcircuit::QMapper::check_connectivity(list<int>& g_front, list<int>& g_pair, list<int>& g_sat, list<int>& g_vio, vector<int>& frozen, Circuit& dgraph)
{
    //(1) Single-qubit gate
    vector<int> g_front_erase;
    for(auto& gateid : g_front)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
            //1. cnot gate
        if(control != -1) continue; 
            //2. single-qubit gate
        else{
            g_front_erase.push_back(gateid);
            g_sat.push_back(gateid);
        }
    }
    for(auto gateid : g_front_erase)
        g_front.remove(gateid);
        
    //(2) CNOT ------
    bool complete_g_pair = false;
    vector<int> g_pair_erase;
    for(auto& gateid : g_pair)
    {
        int control   = dgraph.nodeset[gateid].control;
        int target    = dgraph.nodeset[gateid].target;
        int p_control = layout.log2phy(control);
        int p_target  = layout.log2phy(target );
        //0. gateid=-1 -> skip
        if(gateid == -1)
            g_pair_erase.push_back(gateid);
        else
        {
            //1. CNOT
            if(coupling_graph.dist[p_control][p_target] == 1)
            {
                //erase
                g_pair_erase.push_back(gateid);
                g_sat.push_back(gateid);

                //NOTE: added by dhkim, 2022.04.27
                //FIXME: this is just a quick-fix, might not resolve the fundamental issue
                //- issue: one gate was added to both g_sat[] and g_vio[]
                if (find(g_vio.begin(), g_vio.end(), gateid) != g_vio.end()) {
                    // remove gate from g_vio if the gate is executable
                    g_vio.remove(gateid);
                }

                //push to Final circuit
                complete_g_pair = true;
            }
            else //g_vio
                g_vio.push_back(gateid);
        }
        
    }
    for(auto gateid : g_pair_erase)
        g_pair.remove(gateid);

    //erase duplication
    g_sat.sort();
    g_sat.unique();
    g_vio.sort();
    g_vio.unique();

    return complete_g_pair;
}

void Qcircuit::QMapper::update_bridge_gates(list<int>& g_bridge, list<int>& g_vio, vector<int>& frozen_constraints, Circuit& dgraph)
{
    g_bridge.clear();
    for(auto& gateid : g_vio)
    {
        int control   = dgraph.nodeset[gateid].control;
        int target    = dgraph.nodeset[gateid].target;
        int p_control = layout.log2phy(control);
        int p_target  = layout.log2phy(target );
        if(frozen_constraints[p_control] == 1 || frozen_constraints[p_target] == 1) continue;
        if(coupling_graph.dist[p_control][p_target] == 2)
        {
            int qb = -1;
            for(int i=0; i<positions; i++)
            {
                if(coupling_graph.dist[p_control][i] == 1)
                {
                    if(coupling_graph.dist[p_target][i] == 1)
                    {
                        if(frozen_constraints[i] == 1) continue;
                        if(layout.phy2log(i) < nqubits) //not contain empty qubit
                            qb = i;
                    }
                }
            }
            if(qb != -1)
                g_bridge.push_back(gateid);
        }
    }
}

void Qcircuit::QMapper::generate_swap_candi_list(list<int>& g_vio, vector< pair<pair<int, int>, int> >& swap_candi_list, vector<int>& frozen_constraints, Circuit& dgraph)
{
    //swap candi list
    for(auto& gateid : g_vio)
    {
        int control   = dgraph.nodeset[gateid].control;
        int target    = dgraph.nodeset[gateid].target;
        int p_control = layout.log2phy(control);
        int p_target  = layout.log2phy(target );
        
        for(int i=0; i<coupling_graph.node_size; i++)
        {
            if(layout.phy2log(i) >= nqubits) continue;
            for(auto& p : vector<int>{p_control, p_target} )
            {
                if(coupling_graph.dist[p][i] == 1)
                {
                    if(frozen_constraints[p]==1 || frozen_constraints[i]==1) continue;
                    int SWAP_effect = cal_effect(control, target, i, p);
                    if(SWAP_effect <= 0) continue;
                    
                    int q_phy = layout.phy2log(p);
                    int q_i   = layout.phy2log(i);
                    int Glist_length_phy  = Glist[q_phy].size();
                    int Glist_length_i    = Glist[q_i  ].size();
                    //considering Glist_length
                    if(Glist_length_phy >= Glist_length_i)
                        //phy -> control; i -> target
                        swap_candi_list.push_back(make_pair(make_pair(p, i), gateid));
                    else
                        //i -> control; phy -> target
                        swap_candi_list.push_back(make_pair(make_pair(i, p), gateid));
                }
            }
        }
    }
}

void Qcircuit::QMapper::generate_mov_candi_list(list<int>& g_vio, vector< pair<pair<int, int>, int> >& mov_candi_list, vector<int>& frozen_constraints, Circuit& dgraph)
{
    //mov candi list
    for(auto& gateid : g_vio)
    {
        int control   = dgraph.nodeset[gateid].control;
        int target    = dgraph.nodeset[gateid].target;
        int p_control = layout.log2phy(control);
        int p_target  = layout.log2phy(target );
        
        for(int i=0; i<coupling_graph.node_size; i++)
        {
            if(layout.phy2log(i) < nqubits) continue; //just consider empty qubits
            for(auto& p : vector<int>{p_control, p_target} )
            {
                if(coupling_graph.dist[p][i] == 1)
                {
                    if(frozen_constraints[p]==1 || frozen_constraints[i]==1) continue;
                    int MOV_effect = cal_effect(control, target, i, p);
                    if(MOV_effect <= 0) continue;
                    if(!check_connected_graph(p, i)) continue;

                    mov_candi_list.push_back(make_pair(make_pair(i, p), gateid));
                    
                }
            }
        }
    }
}

bool Qcircuit::QMapper::check_connected_graph(int p_original, int p_new)
{
    int connected = 0;
    for(int qi=0; qi<nqubits; qi++)
    {
        int temp_connected = 0;
        int pi = layout.log2phy(qi);
        if(pi == p_original)
            pi = p_new;
        for(int qj=0; qj<nqubits; qj++)
        {
            if(qi==qj) continue;
            int pj = layout.log2phy(qj);
            if(pj == p_original)
                pj = p_new;
            if(coupling_graph.dist[pi][pj] == 1)
                temp_connected++;
        }
        if(temp_connected != 0) 
            connected++;
    }
    if(connected == nqubits)
        return true;
    else return false;
}

bool compare_priority_gates(pair<int, int> a, pair<int, int> b)
{
    return a.second > b.second;
}

void Qcircuit::QMapper::update_frozen(vector<int>& frozen)
{
    for(int i=0; i<positions; i++)
        frozen[i] = 0;
}

void Qcircuit::QMapper::sort_priority_gates(vector< pair<pair<int, int>, pair<int, int>> >& g_temp, list<int>& g_scheduled, vector<int>& frozen_duration)
{
    g_scheduled.clear();
    vector<vector<pair<int, int>> > temp_gates(2); //0: two-qubit gate, 1: single-qubit gate
    for(auto& gate : g_temp) // gateid, nqubit, kinds, D_length
    {
        if(gate.second.first == 2) // two-qubit gate
            temp_gates[0].push_back(make_pair(gate.first.first, gate.second.second));
        else //single-qubit gate
            temp_gates[1].push_back(make_pair(gate.first.first, gate.second.second));
    }

    for(int i=0; i<2; i++)
    {
        sort(temp_gates[i].begin(), temp_gates[i].end(), compare_priority_gates);
        for(auto& j : temp_gates[i])
            g_scheduled.push_back(j.first);
    }
}
void Qcircuit::QMapper::generate_scheduled_duration(list<int>& gates_list, vector< pair<pair<int, int>, pair<int, int>> >& g_temp, vector<int>& frozen_duration, Circuit& dgraph)
{
    g_temp.clear();

    for(auto& gateid : gates_list)
    {
        const int q_control = dgraph.nodeset[gateid].control;
        const int q_target  = dgraph.nodeset[gateid].target;
        if(q_control == -1) //single-qubit gate
        {
            if(frozen_duration[q_target] == 0) //constraint (@duration) is not exist
            {
                const int Glist_length = Glist[q_target].size();
                g_temp.push_back( make_pair( make_pair(gateid, q_target), make_pair(1, Glist_length) ) ); //(gateid, qubit), (two?or single, length of Glist)
            }
        }
        else //two-qubit gate
        {
            if(frozen_duration[q_control] == 0 && frozen_duration[q_target] == 0) //constraint (@duration) is not exist
            {
                const int Glist_length_c = Glist[q_control].size();
                const int Glist_length_t = Glist[q_target].size();

                int Glist_length_max, q_max;
                int Glist_length_min, q_min;
                if(Glist_length_c >= Glist_length_t)
                {
                    Glist_length_max = Glist_length_c;
                    q_max = q_control;
                    Glist_length_min = Glist_length_t;
                    q_min = q_target;
                }
                else
                {
                    Glist_length_max = Glist_length_t;
                    q_max = q_target;
                    Glist_length_min = Glist_length_c;
                    q_min = q_control;
                }
                g_temp.push_back( make_pair( make_pair(gateid, q_max),  make_pair(2, Glist_length_c + Glist_length_t) ) ); //(gateid, qubit), (two?or single, length of Glist)
            }
        }
    }
}
void Qcircuit::QMapper::update_frozen_duration(vector<int>& frozen_duration)
{
    for(int qi=0; qi<positions; qi++)
    {
        if(frozen_duration[qi] > 0)
            frozen_duration[qi] -= 1;
    }
}

void Qcircuit::QMapper::generate_freq_same(vector< vector<int>>& same_freq_qubits)
{
    for(int i=0; i<positions; i++)
    {
        vector<int> temp_freq;
        const int f_i = coupling_graph.nodeset[i].getweight();
        temp_freq = same_freq_group[f_i];
        temp_freq.erase( remove(temp_freq.begin(), temp_freq.end(), i)  );

        same_freq_qubits.push_back(temp_freq);
    }
}
void Qcircuit::QMapper::generate_freq_relations(vector< pair<vector<int>, vector<int>> >& freq_relations)
{
    for(int i=0; i<positions; i++)
    {
        vector<int> small_freq, large_freq;
        //p = i
        const int f_i = coupling_graph.nodeset[i].getweight();
        for(int j=0; j<positions; j++)
        {
            //p = j
            if(coupling_graph.dist[i][j] == 1)
            {
                const int f_j = coupling_graph.nodeset[j].getweight();
                if(f_i > f_j)
                    small_freq.push_back(j);
                else
                    large_freq.push_back(j);
            }   
        }
        freq_relations.push_back( make_pair( small_freq, large_freq) );
    }
}

void Qcircuit::QMapper::generate_scheduled_frequency(list<int>& gates_list, vector<int>& frozen_frequency, Circuit& dgraph)
{
    list<int> temp_gates_list;
    GATETYPE f1_gatetype, f2_gatetype, f3_gatetype; //U, RX, RZ, CNOT, X, Y, Z, H, S, SDG, T, TDG

    //initiallize
    for(int i=0; i<positions; i++)
        frozen_frequency[i] = -1;
    int i=0;
    for(auto gateid : gates_list)
    {
        const int q_control = dgraph.nodeset[gateid].control;
        const int q_target  = dgraph.nodeset[gateid].target;
        const int p_target  = layout.log2phy(q_target );
        const int f_target  = coupling_graph.nodeset[p_target].getweight();
        const int p_control = layout.log2phy(q_control);
        const int f_control = coupling_graph.nodeset[p_control].getweight();

        // 1. two-qubit gate
        if(q_control != -1) 
        {
            int frozen_control = frozen_frequency[p_control];
            int frozen_target  = frozen_frequency[p_target ];
            if(frozen_control == -1 && frozen_target == -1)
            {
                    // 1-1. f_target (high) -> interact = 1; f_control (low) -> original = 0;
                if(f_control > f_target) 
                {
                    //check neigborhood 
                    int pn_size1  = 0;
                    int pn_check1 = 0;
                    int pn_size2  = 0;
                    int pn_check2 = 0;
                    for(auto& p_i : same_freq_group[f_control]) //-> park (2)
                    {
                        if(coupling_graph.dist[p_target][p_i] == 1 && p_i != p_control)
                        {
                            pn_size1++;
                            if(frozen_frequency[p_i]==-1 || frozen_frequency[p_i]==2)
                                pn_check1++;
                        }
                    }
                    for(auto& p_i : same_freq_group[f_target]) //-> original (0)
                    {
                        if(coupling_graph.dist[p_control][p_i] == 1 && p_i != p_target)
                        {
                            pn_size2++;
                            if(frozen_frequency[p_i]==-1 || frozen_frequency[p_i]==0)
                                pn_check2++;
                        }
                    }
                    //if not violated -> push
                    if(pn_size1 == pn_check1 && pn_size2 == pn_check2)
                    {
                            //gate push
                        temp_gates_list.push_back(gateid);
                        frozen_frequency[p_control] = 0; //original (lower)
                        frozen_frequency[p_target]  = 1; //interactiong (higher->lower)

                        //neighbor frozen_frequency update
                        for(auto& p_i : same_freq_group[f_control])
                            if(coupling_graph.dist[p_target][p_i] == 1 && p_i != p_control)
                                frozen_frequency[p_i] = 2;    //parking (for not interaction)
                        for(auto& p_i : same_freq_group[f_target])
                            if(coupling_graph.dist[p_control][p_i] == 1 && p_i != p_target)
                                frozen_frequency[p_i] = 0;    //original (for not interaction)
                    }
                }
                    // 1-2. f_control (high) -> interact = 1; f_target (low) -> original = 0;
                else                     
                {
                    //check neigborhood 
                    int pn_size1  = 0;
                    int pn_check1 = 0;
                    int pn_size2  = 0;
                    int pn_check2 = 0;
                    for(auto& p_i : same_freq_group[f_target]) //->park (2)
                    {
                        if(coupling_graph.dist[p_control][p_i] == 1 && p_i != p_target)
                        {
                            pn_size1++;
                            if(frozen_frequency[p_i]==-1 || frozen_frequency[p_i]==2)
                                pn_check1++;
                        }
                    }
                    for(auto& p_i : same_freq_group[f_control]) //->original (0)
                    {
                        if(coupling_graph.dist[p_target][p_i] == 1 && p_i != p_control)
                        {
                            pn_size2++;
                            if(frozen_frequency[p_i]==-1 || frozen_frequency[p_i]==0)
                                pn_check2++;
                        }
                    }
                    //if not violated -> push
                    if(pn_size1 == pn_check1 && pn_size2 == pn_check2)
                    {
                            //gate push
                        temp_gates_list.push_back(gateid);
                        frozen_frequency[p_target ] = 0; //original (lower)
                        frozen_frequency[p_control]  = 1; //interactiong (higher->lower)

                        //neighbor frozen_frequency update
                        for(auto& p_i : same_freq_group[f_target])
                            if(coupling_graph.dist[p_control][p_i] == 1 && p_i != p_target)
                                frozen_frequency[p_i] = 2;    //parking (for not interaction)
                        for(auto& p_i : same_freq_group[f_control])
                            if(coupling_graph.dist[p_target][p_i] == 1 && p_i != p_control)
                                frozen_frequency[p_i] = 0;    //original (for not interaction)
                    }

                }
            }
        }
        // 2. single-qubit gate
        else
        {
            int frozen_target = frozen_frequency[p_target ];
            if(frozen_target == -1 || frozen_target == 3)
            {
                int pn_size  = 0;
                int pn_check0 = 0;
                int pn_check3 = 0;
                //check same_freq_group
                for(auto& p_i : same_freq_group[f_target])
                {
                    pn_size++;
                    if(frozen_frequency[p_i]==-1)
                        pn_check0++;
                    if(frozen_frequency[p_i]==3)
                        pn_check3++;
                }
                //if not voilated -> push
                if(pn_size == (pn_check0 + pn_check3))
                {
                    if(pn_check0 == pn_size)
                    {
                        temp_gates_list.push_back(gateid);
                        for(auto& p_i : same_freq_group[f_target])
                            frozen_frequency[p_i] = 3;
                        if(f_target == 1)
                            f1_gatetype = dgraph.nodeset[gateid].type;
                        else if(f_target == 2)
                            f2_gatetype = dgraph.nodeset[gateid].type;
                        else
                            f3_gatetype = dgraph.nodeset[gateid].type;
                    }
                    else
                    {
                        //check gatetype
                        GATETYPE gateid_type = dgraph.nodeset[gateid].type;
                        if(f_target == 1 && (gateid_type == f1_gatetype) )
                            temp_gates_list.push_back(gateid);
                        else if(f_target == 2 && (gateid_type == f2_gatetype) )
                            temp_gates_list.push_back(gateid);
                        else if(f_target == 3 && (gateid_type == f3_gatetype) )
                            temp_gates_list.push_back(gateid);
                    }
                }
            }
        }
        i++;
    }
    gates_list.clear();
    gates_list = temp_gates_list;
}


void Qcircuit::QMapper::update_frozen_frequency(vector<int>& frozen_frequency)
{
    for(int j=0; j<positions; j++)
        frozen_frequency[j] = -1;
}


void Qcircuit::QMapper::print_gates(list<int>& gates, int gate_type) //gate_type) 0: front, 1: pair, 2: satisfied, 3: violated, 4: bridge, 5: scheduled
{
    if(gate_type == 0)
        cout << "- g_front[]: "  << gates << endl;
    else if(gate_type == 1)
        cout << "- g_pair[]: "   << gates << endl;
    else if(gate_type == 2)
        cout << "- g_sat[]: "    << gates << endl;
    else if(gate_type == 3)
        cout << "- g_vio[]: "    << gates << endl;
    else if(gate_type == 4)
    {
        cout << "--------------------------" << endl;
        cout << "- g_bridge[]: " << gates << endl;
    }
    else if(gate_type == 5)
        cout << "- g_scheduled[]: " << gates << endl;
    else
        cout << "* this gate list is not exist" << endl;
}
void Qcircuit::QMapper::print_temp_scheduled_gates(vector<pair<pair<int,int>, pair<int,int>> >& g_temp)
{
    cout << "- g_scheduled (for duration) : \n" << g_temp << endl;
}

void Qcircuit::QMapper::print_gates_spec(list<int>& gates, Circuit& dgraph)
{
    for(auto gateid : gates)
    {
        const int q_control = dgraph.nodeset[gateid].control;
        const int q_target  = dgraph.nodeset[gateid].target;
        const int p_control = layout.log2phy(q_control);
        const int p_target  = layout.log2phy(q_target );
        const int f_control = coupling_graph.nodeset[p_control].getweight();
        const int f_target  = coupling_graph.nodeset[p_target].getweight();
        cout << "* gateid: " << gateid << endl;
        cout << "- control: " << "q" << setw(2) << q_control << ", p" << setw(3) << p_control << ", f" << setw(2) << f_control << endl;
        cout << "- target : " << "q" << setw(2) << q_target  << ", p" << setw(2) << p_target  << ", f" << setw(2) << f_target  << endl;
    }
}

void Qcircuit::QMapper::print_frozen(vector<int>& frozen_vector, int frozen_type)
{
    if(frozen_type == 0)
        cout << "- frozen:                " << frozen_vector << endl;
    else if(frozen_type == 1)
        cout << "- frozen_duration:  " << frozen_vector << endl;
    else if(frozen_type == 2)
        cout << "- frozen_frequency: " << frozen_vector << endl;
    else if(frozen_type == 3)
        cout << "- frozen_constraints:    " << frozen_vector << endl;
    else
        cout << "* this frozen is not exist" << endl;
}

void Qcircuit::QMapper::print_Glist()
{
    cout << "\n* Glist ** " << endl;
    for(int n=0; n<Glist.size(); n++)
    {
        cout << "q" << setw(2) << n << "(p" << setw(2) << layout.log2phy(n) << "): size: " << setw(4) << Glist[n].size() << " ";
        std::list<pair<int, int> >::const_iterator it = Glist[n].begin();
        int num = 13;
        while(num !=0 && it != Glist[n].end())
        {
            if(it->first != -1)
                cout << "-> (" << setw(3) << it->first << ", " << it->second << ")";
            else
                cout << "-> - - - - ";
            it++;
            num--;
        }
        cout << endl;
    }
    cout << "* log->phy" << endl;
    layout.print_log2phy();
    //cout << "* phy->log" << endl;
    //layout.print_phy2log();
}

void Qcircuit::QMapper::print_candi_list(vector< pair< pair<int,int>, int> >& candi_list, int candi_type)
{
    if(candi_type == 0)
    {
        cout << "- swap candi list : ----------------" << endl;
        cout << candi_list << endl;
    }
    else if(candi_type == 1)
    {
        cout << "- mov candi list : -----------------" << endl;
        cout << candi_list << endl;
    }
    else
        cout << "* this candi_list is not exist" << endl;
}

void Qcircuit::QMapper::print_addGates(int gateid, int q_control, int q_target, int p_control, int p_target)
{
    cout << "* gateid: " << gateid << endl;
    cout << "- control: " << "q" << setw(2) << q_control << ", p" << setw(2) << p_control << endl; //<< ", f" << setw(2) << f_control << endl;
    cout << "- target : " << "q" << setw(2) << q_target  << ", p" << setw(2) << p_target  << endl; //<< ", f" << setw(2) << f_target  << endl;

}

void Qcircuit::QMapper::print_candi_pair(pair<int,int>& Pair, int gateid, int max_cost, int candi_type) //candi_type) 0: swap, 1: mov
{
    if(candi_type == 0)
    {
        cout << "SWAP: "          << Pair     << endl;
        cout << "gateid: "        << gateid   << endl;
        cout << "max_swap_cost: " << max_cost << endl;
    }
    else if(candi_type == 1)
    {
        cout << "MOV: "           << Pair     << endl;
        cout << "gateid: "        << gateid   << endl;
        cout << "max_mov_cost: "  << max_cost << endl;
    }
    else
        cout << "* this candi_type is not exist" << endl;
}
