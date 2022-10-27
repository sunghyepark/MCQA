// mappingFunciton.cpp

#include "circuit.h"

using namespace std;
using namespace Qcircuit;


bool compare_by_cost(pair< pair<int,int>, pair<int,int> > p1, pair< pair<int,int>, pair<int,int> > p2)
{
    return p1.second.second > p2.second.second; //original
}

int Qcircuit::QMapper::cal_effect(const int q1, const int q2, const int p1, const int p2)
{
    int p_control = layout.log2phy(q1);
    int p_target  = layout.log2phy(q2);
    int dist_before = coupling_graph.dist[p_control][p_target];

    if(p_control == p1 || p_control == p2)
        p_control = (p_control == p1) ? p2 : p1;
    if(p_target  == p1 || p_target  == p2)
        p_target  = (p_target  == p1) ? p2 : p1;

    int dist_after = coupling_graph.dist[p_control][p_target];
    int effect = dist_before - dist_after;

    return effect;
}

double Qcircuit::QMapper::cal_MCPE(const pair<int, int> p, Circuit& dgraph, bool isswap)
{
    int q1 = layout.phy2log(p.first );
    int q2 = layout.phy2log(p.second);

    double MCPE = 0;
    int dist = 0;
    double power = 1.0;

    //cout << "** cal_MCPE [ p" << p.first << "(q" << q1 << ") ,p" << p.second << "(q" << q2 << ") ] **" << endl;

    //for mov
    if(!isswap)
        if(q1 > q2)
            q1 = q2;

    //for q1
    if(q1 < nqubits)
    {
        for(auto& node : Glist[q1])
        {
            int gateid = node.first;
            int control, target;
            int effect;
            if(gateid != -1)
            {
                Gate &g = dgraph.nodeset[gateid];
                control = g.control;
                target  = g.target;
                if(control != -1)
                {
                    effect = cal_effect(control, target, p.first, p.second);
                    if(effect < 0) break;
                    if(power < 0.00000001) break;
                    MCPE += effect * power;
                    power*= param_alpha;
                    dist++;
                }
                else
                {
                    power *= param_alpha;
                    dist++;
                }
            }
        }
    }

    //for swap
    if(isswap)
    {
        //for q2
        dist = 0;
        power = 1.0;
        if(q2 < nqubits)
        {
            for(auto& node : Glist[q2])
            {
                int gateid = node.first;
                int control, target;
                int effect;
                if(gateid != -1)
                {
                    Gate &g = dgraph.nodeset[gateid];
                    control = g.control;
                    target  = g.target;
                    if(control != -1)
                    {
                        effect = cal_effect(control, target, p.first, p.second);
                        if(effect < 0) break;
                        if(power < 0.00000001) break;
                        MCPE += effect * power;
                        power*= param_alpha;
                        dist++;
                    }
                    else
                    {
                        power *= param_alpha;
                        dist++;
                    }
                }
            }
        }
    }
    //cout << "-> MCPE: " << MCPE << endl << endl;
    return MCPE*10;
}

void Qcircuit::QMapper::find_max_cost(pair<int, int>& SWAP, int& gateid, double& max_swap_cost,
                                            vector< pair< pair<int,int> , pair<int,double> > >& v, list<int>& non_exec_gates,
                                            vector< pair< pair<int,int> , pair<int,double> > >& history)
{
    sort(v.begin(), v.end(), compare_by_cost);
    int index = 0;
    pair< pair<int, int>, pair<int, double> >& ref = v[index];
    for(auto pair : history)
    {
        if(pair == ref)
            index++;
    }

    if(index >= v.size()) index = v.size()-1;
    //cout << "sorted" << endl;
    //cout << v << endl;

    SWAP          = v[index].first;
    gateid        = v[index].second.first;
    max_swap_cost = v[index].second.second;
}

void Qcircuit::QMapper::add_cnot(int c_qubit, int t_qubit, Circuit& graph)
{
#if GATESET_IIC_JKU || GATESET_DAC2022
    Qcircuit::Gate cnot;
    
    cnot.id = node_id;
    cnot.control = c_qubit;
    cnot.target = t_qubit;
    cnot.type = GATETYPE::CNOT;
    
    graph.nodeset.push_back(cnot);

    add_cnot_num++;
    node_id++;
#elif GATESET_QUTECH
    Qcircuit::Gate rym90, ry90;
    
    rym90.id = node_id;
    rym90.target = t_qubit;
    rym90.control = -1;
    rym90.type = GATETYPE::RY;
    snprintf(rym90.output_type, 127, "(%.0f)", -90.);
    graph.nodeset.push_back(rym90);
    node_id++;

    add_cz(c_qubit, t_qubit, graph);

    ry90.id = node_id;
    ry90.target = t_qubit;
    ry90.control = -1;
    ry90.type = GATETYPE::RY;
    snprintf(ry90.output_type, 127, "(%.0f)", 90.);
    graph.nodeset.push_back(ry90);
    node_id++;
#endif
}

#if GATESET_QUTECH
void Qcircuit::QMapper::add_cz(int c_qubit, int t_qubit, Circuit& graph){
    Qcircuit::Gate cz;
    
    cz.id = node_id;
    cz.control = c_qubit;
    cz.target = t_qubit;
    cz.type = GATETYPE::CZ;
    snprintf(cz.output_type, 127, "CZ");
    graph.nodeset.push_back(cz);
    add_cz_num++;
    node_id++;
}
#endif

void Qcircuit::QMapper::add_swap(int q1, int q2, Circuit& graph)
{
    add_cnot(q1, q2, graph);
    add_cnot(q2, q1, graph);
    add_cnot(q1, q2, graph);

    add_swap_num++;
}

void Qcircuit::QMapper::add_mov(int q1, int q2, Circuit& graph)
{
    add_cnot(q1, q2, graph);
    add_cnot(q2, q1, graph);


    add_mov_num++;
}

void Qcircuit::QMapper::add_bridge(int ps, int pt, vector<int>& frozen_constraints, Circuit& graph)
{
    int pb;
    int qs, qt, q;
    /*ps: control, pt: target*/
    for(int i=0; i<positions; i++)
    {
        if(coupling_graph.dist[ps][i] == 1)
        {
            if(coupling_graph.dist[pt][i] == 1)
            {
                if(layout.phy2log(i) >= nqubits) continue;
                if(frozen_constraints[i] == 1) continue;
                pb = i;
            }
            else continue;
        }
        else continue;
    }

    qs = layout.phy2log(ps);
    qt = layout.phy2log(pt);
    q  = layout.phy2log(pb);
    
    add_cnot(q,  qt, graph);
    add_cnot(qs, q, graph);
    add_cnot(q,  qt, graph);
    add_cnot(qs, q, graph);
    
#if GATESET_IIC_JKU || GATESET_DAC2022
    add_cnot_num--;
#elif GATESET_QUTECH
    add_cz_num--;
#endif
    add_bridge_num++;
}

void Qcircuit::QMapper::add_gate_to_Glist(int gate_mode, int update_nodeid, int q_control, int q_target, int q_middle)
{
    //BRIDGE
    if(gate_mode == 0)
    {
        int bridge_nodeid = update_nodeid;

        int front_control = Glist[q_control].front().first;
        int front_target  = Glist[q_target ].front().first;

        _add_gate_to_Glist(q_control, q_middle, bridge_nodeid);
        _add_gate_to_Glist(q_middle,  q_target, bridge_nodeid);
        _add_gate_to_Glist(q_control, q_middle, bridge_nodeid);
        _add_gate_to_Glist(q_middle,  q_target, bridge_nodeid);
    }
    //SWAP
    else if(gate_mode == 1)
    {
        int swap_nodeid = update_nodeid;
        int q_SWAP_first  = q_control;
        int q_SWAP_second = q_target;
        _add_gate_to_Glist(q_SWAP_first,  q_SWAP_second, swap_nodeid);
        _add_gate_to_Glist(q_SWAP_second, q_SWAP_first,  swap_nodeid);
        _add_gate_to_Glist(q_SWAP_first,  q_SWAP_second, swap_nodeid);
    }
    //MOV
    else
    {
        int mov_nodeid = update_nodeid;
        int q_MOV_first  = q_control;
        int q_MOV_second = q_target;
#if GATESET_IIC_JKU || GATESET_DAC2022
        _add_gate_to_Glist(q_MOV_first,  q_MOV_second, mov_nodeid);
        _add_gate_to_Glist(q_MOV_second, q_MOV_first,  mov_nodeid);
#elif GATESET_QUTECH
        _add_gate_to_Glist(q_MOV_second, q_MOV_first,  mov_nodeid);
        _add_gate_to_Glist(q_MOV_first,  q_MOV_second, mov_nodeid);
#endif
    }
}

#if GATESET_IIC_JKU || GATESET_DAC2022
void Qcircuit::QMapper::_add_gate_to_Glist(int q1, int q2, int &gateid)
{
    auto it1 = Glist[q1].begin(); 
    auto it2 = Glist[q2].begin(); 

    if(it1 != Glist[q1].end())
        if(Glist[q1].front().first == -1)
            it1++;
    if(it2 != Glist[q2].end())
        if(Glist[q2].front().first == -1)
            it2++;

    gateid--;
    Glist[q2].insert( it2, (make_pair(gateid, 2)) ); 
    Glist[q2].insert( it2, (make_pair(-1, 1))     ); 
    Glist[q1].insert( it1, (make_pair(gateid, 2)) ); 
    Glist[q1].insert( it1, (make_pair(-1, 1))     ); 
}
#elif GATESET_QUTECH
void Qcircuit::QMapper::_add_gate_to_Glist(int q1, int q2, int &gateid)
{
    auto it1 = Glist[q1].begin(); 
    auto it2 = Glist[q2].begin(); 

    if(it1 != Glist[q1].end())
        if(Glist[q1].front().first == -1)
            it1++;
    if(it2 != Glist[q2].end())
        if(Glist[q2].front().first == -1)
            it2++; 
    
    gateid--;
    Glist[q2].insert( it2, (make_pair(gateid-2, 1)) ); 
    gateid--;
    Glist[q2].insert( it2, (make_pair(gateid, 2)) ); 
    Glist[q2].insert( it2, (make_pair(-1, 1))     ); 
    Glist[q1].insert( it1, (make_pair(gateid, 2)) ); 
    Glist[q1].insert( it1, (make_pair(-1, 1))     ); 
    gateid--;
    Glist[q2].insert( it2, (make_pair(gateid+2, 1)) ); 
}
#endif
