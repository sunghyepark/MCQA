// initial_mapping.cpp

#include "circuit.h"
#include <random>

using namespace std;
using namespace Qcircuit;

bool compare_interaction_weight(pair<int, int> a, pair<int, int> b)
{
    return a.second > b.second;
}

void Qcircuit::QMapper::initial_mapping(ARCHITECTURE archi, bool prepro, bool BRIDGE_MODE, INITIALTYPE initial_type)
{
    //*** Initial mapping generate
    if(!prepro)
    {
        for(int i=0; i<positions; i++)
        {
            layout_L[i] = i; //logical -> physical
            qubit_Q[i]  = i; //physical -> logical
        }
        return;
    }
    
    //*** Coupling graph generate
    select_coupling_graph(archi);
    coupling_graph.build_dist_table();

    //Initial mapping mode?
    int q_num;
    switch(initial_type)
    {
        case INITIALTYPE::IDENTITY_MAPPING:
            cout << "Initial type: identity mapping" << endl;
            identical_mapping();
            break;
        case INITIALTYPE::RANDOM_MAPPING:
            cout << "Initial type: random mapping" << endl;
            random_mapping();
            break;
        case INITIALTYPE::GRAPH_MATCHING_PROPOSED:
            //cout << "Initial type: graph matching (proposed)" << endl;
            graph_matching_processing_proposed(archi);
            q_num = nqubits;
            for(int i=0; i<positions; i++)
            {
                if(qubit_Q[i] == -1)
                {
                    qubit_Q[i] = q_num;
                    layout_L[q_num] = i;
                    q_num++;
                }
            }
            break;
    }

    // to save the result of initial mapping
    initial_layout = Layout(layout_L);
    // to use in main mapping process
    layout = Layout(layout_L);

    // # Final initial mapping 
    /*
    cout << "* logical -> physical" << endl;
    layout.print_log2phy();
    cout << "* physical -> logical" << endl;
    layout.print_phy2log();
    */
}

void Qcircuit::QMapper::graph_matching_processing_proposed(ARCHITECTURE archi)
{
    // (0) initiallize ========================================================
    //initiallize layout -> -1
    for(int i=0; i<positions; i++)
        layout_L[i] = -1; //logical -> physical
    for(int i=0; i<positions; i++)
        qubit_Q[i]  = -1; //physical -> logical

    //CNOT set
    make_CNOT(true);
    
    //paramter setting
    int nodesetSize_cnot = Dgraph_cnot.nodeset.size();
    int nodesetSize      = Dgraph.nodeset.size();
    int param = 9;
    int n = nodesetSize_cnot / param;
    if(nqubits <= 10 || nqubits == 20)
        n = 1;
    if( (nqubits == 15 && nodesetSize <= 5000) || nqubits == 11)
        n = param;
    
    // make interaction graph
    interactionGraph = generate_interaction_graph(true, n); //G_i
    interactionGraph.build_dist_table();

    // (1) grpah matching start =================================================
    int qc = interactionGraph.generate_graph_center();
    int pc = coupling_graph.generate_graph_center();
    qubit_Q[pc ] = qc;
    layout_L[qc] = pc; 

    queue<int> igraph_qc_queue;
        generate_bfs_queue(igraph_qc_queue, interactionGraph, qc, false);
        
        while(!igraph_qc_queue.empty())
        {
            int current_qc = igraph_qc_queue.front();

            //(1-1) qc mapping ------------------------------------------------------
            //check layout
            if(layout_L[current_qc] == -1)
            {
                vector<int> ref_loc, candi_loc;
                    //make ref_loc
                make_ref_loc(ref_loc, interactionGraph, current_qc, false);
                    //make candi_loc
                make_candi_loc(current_qc, candi_loc, ref_loc, coupling_graph, interactionGraph, false);
                sort_degree(candi_loc, coupling_graph);
                pc = candi_loc[0];
                    //map current_qc <-> pc
                qubit_Q[pc         ] = current_qc;
                layout_L[current_qc] = pc;
            }
            
            //(1-2) dist_1_mapping --------------------------------------------------
            //make dist_1_group
            queue<int> current_qc_queue;
            vector<int> dist_1_group;
            generate_bfs_queue(current_qc_queue, interactionGraph, current_qc, false);
        
        while(!current_qc_queue.empty())
        {
            int front = current_qc_queue.front();
            if(current_qc != front)
                if(interactionGraph.dist[current_qc][front] == 1)
                    dist_1_group.push_back(front);
            current_qc_queue.pop();
        }
        
        //make candidate_loc
        vector<int> candi_loc, ref_loc;
        ref_loc.push_back(layout_L[current_qc]);
        make_candi_loc(current_qc, candi_loc, ref_loc, coupling_graph, interactionGraph, false);
        sort_degree(candi_loc, coupling_graph);

        //dist_1_mapping
        for(int idx = 0; idx < dist_1_group.size(); idx++)
        {
            int current_dist_1 = dist_1_group[idx];
            if(layout_L[current_dist_1] != -1) continue;
            if(idx == 0)
            {
                int first_candiloc = candi_loc[0];
                qubit_Q[first_candiloc ] = current_dist_1;
                layout_L[current_dist_1] = first_candiloc;
                candi_loc.erase( remove(candi_loc.begin(), candi_loc.end(), first_candiloc) );
            }
            else
            {
                if(candi_loc.empty()) continue; 
                int candi_layout;
                vector<pair<int,int>> weight_relations;
                for(int j=0; j<idx; j++)
                {
                    //interaction exist
                    if(interactionGraph.dist[current_dist_1][dist_1_group[j]] == 1)
                    {
                        int w_interaction = get_interaction_weight(current_dist_1, dist_1_group[j]);
                        weight_relations.push_back( make_pair(dist_1_group[j], w_interaction) );
                    }
                }

                //1. no interaction w/ other dist1 nodes
                if(weight_relations.empty())
                {
                    vector<int> dist_1_layout;
                    for(int i = 0; i < idx; i++)
                        dist_1_layout.push_back( layout_L[dist_1_group[i]] );
                    candi_layout = generate_candi_layout( candi_loc, dist_1_layout, true, true);
                }
                //2. interaction exists w/ other dist1 nodes
                else
                {
                    int w_interaction = get_interaction_weight(current_dist_1, current_qc);
                    weight_relations.push_back( make_pair(current_qc, w_interaction) );
                    sort(weight_relations.begin(), weight_relations.end(), compare_interaction_weight);
                    if(weight_relations[0].first == current_qc)
                    {
                        if(weight_relations.size() > 1)
                        {
                            vector<int> candi_loc2, ref_loc2;
                            int q_first_relations = weight_relations[1].first;
                            ref_loc2.push_back( layout_L[q_first_relations] );
                            make_candi_loc(q_first_relations, candi_loc2, ref_loc2, coupling_graph, interactionGraph, false);
                            candi_layout = generate_candi_layout( candi_loc, ref_loc2, false, true);
                        }
                        else
                        {
                            sort_degree(candi_loc, coupling_graph);
                            candi_layout = candi_loc[0];
                        }
                    }
                    else
                    {
                        vector<int> candi_loc2, ref_loc2;
                        int q_first_relations = weight_relations[0].first;
                        ref_loc2.push_back( layout_L[q_first_relations] );
                        make_candi_loc(q_first_relations, candi_loc2, ref_loc2, coupling_graph, interactionGraph, false);
                        vector<int> _candi_loc;
                        _candi_loc.push_back(pc);
                        candi_layout = generate_candi_layout( candi_loc2, _candi_loc, false, true);
                    }
                }
                
                //qc -> candi_layout
                qubit_Q[candi_layout   ] = current_dist_1;
                layout_L[current_dist_1] = candi_layout;
                auto it = std::find( candi_loc.begin(), candi_loc.end(), candi_layout );
                if(it != candi_loc.end() )
                    candi_loc.erase( remove(candi_loc.begin(), candi_loc.end(), candi_layout) );
            }
        }
        //(1-3) prepare next step -----------------------------------------------
        igraph_qc_queue.pop();
    }

    //(2) frequency matching start ===================================================
    //print_layout(layout_L); //prev layout
    
    //(2-0) layout -> make layoutGraph --------------------------------
    //make layoutGraph
    layoutGraph = make_layoutGraph(true, false);
    layoutGraph.build_dist_table();
    
    //make p_num -> (x, y)
    if(archi == ARCHITECTURE::Surface_17)
    {
        build_coordinate_Surface_17(); //p_to_xy; (pi -> (x,y))
        freq_mod.push_back(3); //idx:0
        freq_mod.push_back(2); //idx:1
        freq_mod.push_back(1); //idx:2
        freq_mod.push_back(2); //idx:3
    }
   
    int score_origin = return_layout_score(layoutGraph, false);
    int score_max = score_origin;
    int best_isReflected, best_move_x, best_move_y;
    for(int move_x = -1; move_x<=1; move_x++)
    {
        for(int move_y = -1; move_y<=1; move_y++)
        {
            for(int isReflected = 0; isReflected<=1; isReflected++)
            {
                if(move_x == 0 && move_y == 0 && isReflected == 0) continue;

                int generated = make_layout(move_x, move_y, isReflected);
                if(generated)
                {
                    new_layoutGraph = make_layoutGraph(true, true);
                    new_layoutGraph.build_dist_table();
                    int score_new = return_layout_score(new_layoutGraph, true);
                    if(score_new > score_max)
                    {
                        score_max = score_new;
                        best_isReflected  = isReflected;
                        best_move_x = move_x;
                        best_move_y = move_y;
                    }
                }
            }
        }
    }
    //(2-2) find best layout --------------------------------------
    //make edge_available_cnot_num
    if(score_max > score_origin)
    {
        int generated = make_layout(best_move_x, best_move_y, best_isReflected);
        for(int i=0; i<positions; i++)
            layout_L[i] = new_layout_L[i]; //logical -> physical
        for(int i=0; i<positions; i++)
            qubit_Q[i ] = new_qubit_Q[i ]; //physical -> logical
    }
}

void Qcircuit::QMapper::identical_mapping()
{
    for(int i=0; i<positions; i++)
    {
        qubit_Q[i]  = i;
        layout_L[i] = i;
    }
}

void Qcircuit::QMapper::random_mapping()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0,positions);
    
    int rand_num;
    int i=0;

    //initiallize layout to -1
    for(int i=0; i<nqubits; i++)
        layout_L[i] = -1; //logical -> physical
    for(int i=0; i<positions; i++)
        qubit_Q[i] = -1;  //physical -> logical

    while(i<positions)
    {
        rand_num = dis(gen);
        if(qubit_Q[rand_num] != -1) continue;
        else
        {
            qubit_Q[rand_num] = i;
            layout_L[i] = rand_num;
            cout << "qubit: " << i << ", layout: " << rand_num << endl;
            i++;
        }
    }
}

int Qcircuit::QMapper::get_interaction_weight(int sourceid, int targetid)
{
    int weight;
    for(auto& e : interactionGraph.edgeset)
    {
        if(e.second.getsourceid() == sourceid && e.second.gettargetid() == targetid
        || e.second.getsourceid() == targetid && e.second.gettargetid() == sourceid)
            weight = e.second.getweight();
    }
    return weight;
}

int Qcircuit::QMapper::generate_candi_layout(vector<int>& candi_loc1, vector<int>& candi_loc2, bool dist_max, bool forward)
{
    int max = -1;
    int min = 10000;
    vector<int> temp_candi_layout;
    
    //(1)
    if(candi_loc1.size() >= candi_loc2.size())
    {
        for(auto Q_candi : candi_loc1)
        {
            //(1-1) calculate temp_dist
            int temp_dist = 0;
            for(int i=0; i<candi_loc2.size(); i++)
            {
                int Q_current    = candi_loc2[i];
                int dist_current = coupling_graph.dist[Q_candi][Q_current];
                temp_dist += dist_current;
            }

            //(1-2) make temp_candi_layout
            if(dist_max == true)
            {
                if(temp_dist > max)
                {
                    max = temp_dist;
                    temp_candi_layout.clear();
                    temp_candi_layout.push_back(Q_candi);
                }
                else if(temp_dist == max)
                    temp_candi_layout.push_back(Q_candi);
            }
            else
            {
                if(temp_dist < min)
                {
                    min = temp_dist;
                    temp_candi_layout.clear();
                    temp_candi_layout.push_back(Q_candi);
                }
                else if(temp_dist == min)
                    temp_candi_layout.push_back(Q_candi);
            }
        }
    }
    //(2)
    else
    {
        for(auto Q_candi : candi_loc1)
        {
            //(2-1) calculate temp_dist
            int temp_dist = 0;
            for(int i=0; i<candi_loc1.size(); i++)
            {
                int Q_current    = candi_loc2[i];
                int dist_current = coupling_graph.dist[Q_candi][Q_current];
                temp_dist += dist_current;
            }
            
            //(2-2) make temp_candi_layout
            if(dist_max == true)
            {
                if(temp_dist > max)
                {
                    max = temp_dist;
                    temp_candi_layout.clear();
                    temp_candi_layout.push_back(Q_candi);
                }
                else if(temp_dist == max)
                    temp_candi_layout.push_back(Q_candi);
            }
            else
            {
                if(temp_dist < min)
                {
                    min = temp_dist;
                    temp_candi_layout.clear();
                    temp_candi_layout.push_back(Q_candi);
                }
                else if(temp_dist == min)
                    temp_candi_layout.push_back(Q_candi);
            }
        }
    }

    //(3) sort temp_candi_layout
    sort_degree(temp_candi_layout, coupling_graph);
    
    //(4) return candi_layout
    if(!forward)
        reverse( temp_candi_layout.begin(), temp_candi_layout.end() );
    return temp_candi_layout[0];
}

void Qcircuit::QMapper::print_layout(const map<int, int> &layout)
{
    cout.setf(ios::left);
    for(auto& it : layout)
        cout << " q" << setw(2) << it.first << "\t->    Q" << setw(2) << it.second << endl;
}

void Qcircuit::QMapper::print_qubit(const map<int, int> &qubit)
{
    cout.setf(ios::left);
    for(auto& it : qubit)
        cout << " Q" << setw(2) << it.first << "\t<-    q" << setw(2) << it.second << endl;
}

void Qcircuit::QMapper::print_layout()
{
    const map<int, int> &layout = layout_L;
    cout.setf(ios::left);
    for(auto& it : layout)
        cout << " q" << setw(2) << it.first << "\t->    Q" << setw(2) << it.second << endl;

}

void Qcircuit::QMapper::print_qubit()
{
    const map<int, int> &qubit = qubit_Q;
    cout.setf(ios::left);
    for(auto& it : qubit)
        cout << " q" << setw(2) << it.first << "\t<-    q" << setw(2) << it.second << endl;

}

int Qcircuit::QMapper::return_freq_idx(int x, int y)
{
    int sum = x+y;
    int mod = sum % 4;
    if(mod < 0) mod+= 4;
    return mod;
}

void Qcircuit::QMapper::return_noninteraction_edges(vector<pair<int, int>>& noninteraction_edges, int Qi, int Qj)
{
    int xi, yi, xj, yj;
    xi = p_to_xy[Qi].first ;
    xj = p_to_xy[Qj].first ;
    yi = p_to_xy[Qi].second;
    yj = p_to_xy[Qj].second;
    int fi_idx = return_freq_idx(xi, yi);
    int fj_idx = return_freq_idx(xj, yj);
    int fi = freq_mod[fi_idx];
    int fj = freq_mod[fj_idx];
    
    //neigbor coordinate
    int Qi_u   = xy_to_p(xi  , yi+1); //i
    int Qi_d   = xy_to_p(xi  , yi-1);
    int Qi_r   = xy_to_p(xi+1, yi  );
    int Qi_l   = xy_to_p(xi-1, yi  );
    int Qi_ul  = xy_to_p(xi-1, yi+1);
    int Qi_ur  = xy_to_p(xi+1, yi+1);
    int Qi_dl  = xy_to_p(xi-1, yi-1);
    int Qi_dr  = xy_to_p(xi+1, yi-1);
    int Qi_2u  = xy_to_p(xi  , yi+2);
    int Qi_2d  = xy_to_p(xi  , yi-2);
    int Qi_2r  = xy_to_p(xi+2, yi  );
    int Qi_2l  = xy_to_p(xi-2, yi  );
    int Qj_u   = xy_to_p(xj  , yj+1); //j
    int Qj_d   = xy_to_p(xj  , yj-1);
    int Qj_r   = xy_to_p(xj+1, yj  );
    int Qj_l   = xy_to_p(xj-1, yj  );
    int Qj_ul  = xy_to_p(xj-1, yj+1);
    int Qj_ur  = xy_to_p(xj+1, yj+1);
    int Qj_dl  = xy_to_p(xj-1, yj-1);
    int Qj_dr  = xy_to_p(xj+1, yj-1);
    int Qj_2u  = xy_to_p(xj  , yj+2);
    int Qj_2d  = xy_to_p(xj  , yj-2);
    int Qj_2r  = xy_to_p(xj+2, yj  );
    int Qj_2l  = xy_to_p(xj-2, yj  );

    if(fi == 1 && fj== 2)
    {
        //1
        if(fj_idx == 1) // fi=1[2] <-> fj=2[1]
        {
            if(Qi_u != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_u)  );
            if(Qi_r != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_r)  );
            if(Qj_l != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_l)  );
            if(Qj_d != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_d)  );
            if(Qi_u != -1 && Qi_2u != -1)
                noninteraction_edges.push_back( make_pair(Qi_u, Qi_2u)  );
            if(Qi_u != -1 && Qi_ul != -1)
                noninteraction_edges.push_back( make_pair(Qi_u, Qi_ul)  );
            if(Qi_r != -1 && Qi_ur != -1)
                noninteraction_edges.push_back( make_pair(Qi_r, Qi_ur)  );
            if(Qi_r != -1 && Qi_2r != -1)
                noninteraction_edges.push_back( make_pair(Qi_r, Qi_2r)  );
            //1-1
            if( abs(yi-yj)==1  )
            {
                if(Qi_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_l)  );
                if(Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_r)  );
                if(Qi_ul != -1 && Qi_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_ul, Qi_l)  );
                if(Qi_l != -1 && Qi_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_2l)  );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qj_l)  );
                if(Qi_u != -1 && Qi_ur != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_ur)  );
                if(Qj_r != -1 && Qj_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_2r)  );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r)  );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr)  );
            }
            //1-2
            else
            {
                if(Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_d)  );
                if(Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_u)  );
                if(Qj_u != -1 && Qj_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_2u)  );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qj_u)  );
                if(Qj_u != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_ul)  );
                if(Qj_d != -1 && Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qi_d)  );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr)  );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_dr)  );
                if(Qi_d != -1 && Qi_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_2d)  );
            } 
        }
        //2
        else            // fi=1[2] <-> fj=2[3]
        {
            if(Qi_l != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_l)  );
            if(Qi_d != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_d)  );
            if(Qj_r != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_r)  );
            if(Qj_u != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_u)  );
            if(Qi_l != -1 && Qi_2l != -1)
                noninteraction_edges.push_back( make_pair(Qi_l, Qi_2l)  );
            if(Qi_l != -1 && Qi_dl != -1)
                noninteraction_edges.push_back( make_pair(Qi_l, Qi_dl)  );
            if(Qi_d != -1 && Qi_dl != -1)
                noninteraction_edges.push_back( make_pair(Qi_d, Qi_dl)  );
            if(Qi_d != -1 && Qi_2d != -1)
                noninteraction_edges.push_back( make_pair(Qi_d, Qi_2d)  );
            //2-1
            if( abs(yi-yj)==1  )
            {
                if(Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_r)  );
                if(Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_l)  );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_dr)  );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr)  );
                if(Qi_r != -1 && Qi_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_2r)  );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qj_l)  );
                if(Qj_l != -1 && Qj_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_2l)  );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r)  );
            }
            //2-2
            else
            {
                if(Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_u) );
                if(Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_d) );
                if(Qi_u != -1 && Qi_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_2u) );
                if(Qi_u != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_ul) );
                if(Qi_l != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_ul) );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qj_u) );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d) );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_dr) );
                if(Qj_d != -1 && Qj_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_2d) );
            } 
        }
    }
    else if(fi == 2 && fj == 1)
    {
        //1
        if(fi_idx == 1) // fi=2[1] <-> fj=1[2]
        {
            if(Qj_u != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_u)  );
            if(Qj_r != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_r)  );
            if(Qi_l != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_l)  );
            if(Qi_d != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_d)  );
            if(Qj_u != -1 && Qj_2u != -1)
                noninteraction_edges.push_back( make_pair(Qj_u, Qj_2u)  );
            if(Qj_u != -1 && Qj_ur != -1)
                noninteraction_edges.push_back( make_pair(Qj_u, Qj_ur)  );
            if(Qj_r != -1 && Qj_ur != -1)
                noninteraction_edges.push_back( make_pair(Qj_r, Qj_ur)  );
            if(Qj_r != -1 && Qj_2r != -1)
                noninteraction_edges.push_back( make_pair(Qj_r, Qj_2r)  );
            //1-1
            if( abs(yi-yj)==1  )
            {
                if(Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_l)  );
                if(Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_r)  );
                if(Qj_u != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_ul)  );
                if(Qj_ul != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj_ul, Qj_l)  );
                if(Qj_l != -1 && Qj_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_2l)  );
                if(Qj_l != -1 && Qi_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qi_l)  );
                if(Qi_r != -1 && Qi_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_2r)  );
                if(Qj_r != -1 && Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qi_r)  );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr)  );
            }
            //1-2
            else
            {
                if(Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_d)  );
                if(Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_u)  );
                if(Qi_u != -1 && Qi_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_2u)  );
                if(Qj_u != -1 && Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qi_u)  );
                if(Qi_u != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_ul)  );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d)  );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr)  );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_dr)  );
                if(Qj_d != -1 && Qj_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_2d)  );
            } 
        }
        //2
        else            // fi=2[3] <-> fj=1[2]
        {
            if(Qj_l != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_l)  );
            if(Qj_d != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_d)  );
            if(Qi_r != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_r)  );
            if(Qi_u != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_u)  );
            if(Qj_l != -1 && Qj_2l != -1)
                noninteraction_edges.push_back( make_pair(Qj_l, Qj_2l)  );
            if(Qj_l != -1 && Qj_dl != -1)
                noninteraction_edges.push_back( make_pair(Qj_l, Qj_dl)  );
            if(Qj_d != -1 && Qj_dl != -1)
                noninteraction_edges.push_back( make_pair(Qj_d, Qj_dl)  );
            if(Qj_d != -1 && Qj_2d != -1)
                noninteraction_edges.push_back( make_pair(Qj_d, Qj_2d)  );
            //2-1
            if( abs(yi-yj)==1  )
            {
                if(Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_r)  );
                if(Qi_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_l)  );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_dr)  );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr)  );
                if(Qj_r != -1 && Qj_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_2r)  );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qj_l)  );
                if(Qi_l != -1 && Qi_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_2l)  );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r)  );
            }
            //2-2
            else
            {
                if(Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_u) );
                if(Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_d) );
                if(Qj_u != -1 && Qj_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_2u) );
                if(Qj_u != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_ul) );
                if(Qj_l != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_ul) );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qj_u) );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d) );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_dr) );
                if(Qi_d != -1 && Qi_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_2d) );
            } 
        }
    }
    else if(fi == 2 && fj == 3)
    {
        //3
        if(fi_idx == 1) // fi=2[1] <-> fj=3[0]
        {
            if(Qi_r != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_r) );
            if(Qj_l != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_l) );
            if(Qi_u != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_u) );
            if(Qj_d != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_d) );
            //3-1
            if( abs(yi-yj)==1  )
            {
                if(Qi_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_l) );
                if(Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_r) );
                if(Qi_l != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_ul) );
                if(Qi_l != -1 && Qi_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_2l) );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l,  Qj_l) );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r) );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr) );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d,  Qj_dr) );
            }
            //3-2
            else
            {
                if(Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_u) );
                if(Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_d) );
                if(Qj_u != -1 && Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qi_u) );
                if(Qj_l != -1 && Qj_dl != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_dl) );
                if(Qj_d != -1 && Qj_dl != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_dl) );
                if(Qj_d != -1 && Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qi_d) );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_dr) );
                if(Qi_d != -1 && Qi_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_2d) );
            } 
        }
        //4
        else            // fi=2[3] <-> fj=3[0]
        {
            if(Qi_l != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_l) );
            if(Qj_r != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_r) );
            if(Qi_d != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_d) );
            if(Qj_u != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_u) );
            //4-1
            if( abs(yi-yj)==1  )
            {
                if(Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_r) );
                if(Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_l) );
                if(Qi_r != -1 && Qi_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_2r) );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr) );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qj_l) );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r) );
                if(Qj_u != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_ul) );
                if(Qj_l != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_ul) );
            }
            //4-2
            else
            {
                if(Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_u) );
                if(Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_d) );
                if(Qi_u != -1 && Qi_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_2u) );
                if(Qi_u != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qi_ul) );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u,   Qj_u) );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d) );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr) );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d,  Qj_dr) );
            } 
        }
    }
    else  //fi == 3 && fj == 2
    {
        //3
        if(fj_idx == 1) // fi=3[0] <-> fj=2[1]
        {
            if(Qj_r != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_r) );
            if(Qj_u != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_u) );
            if(Qi_l != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_l)  );
            if(Qi_d != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_d) );
            //3-1
            if( abs(yi-yj)==1  )
            {
                if(Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_l) );
                if(Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_r) );
                if(Qj_l != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_ul) );
                if(Qj_l != -1 && Qj_2l != -1)
                    noninteraction_edges.push_back( make_pair(Qj_l, Qj_2l) );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l,  Qj_l) );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r) );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr) );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d,  Qi_dr) );
            }
            //3-2
            else
            {
                if(Qi_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_u) );
                if(Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_d) );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u, Qj_u) );
                if(Qi_l != -1 && Qi_dl != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_dl) );
                if(Qi_d != -1 && Qi_dl != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qi_dl) );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d) );
                if(Qj_d != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_dr) );
                if(Qj_d != -1 && Qj_2d != -1)
                    noninteraction_edges.push_back( make_pair(Qj_d, Qj_2d) );
            } 
        }
        //4
        else            // fi=3[0] <-> fj=2[3]
        {
            if(Qj_l != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_l) );
            if(Qi_u != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_u) );
            if(Qj_d != -1)
                noninteraction_edges.push_back( make_pair(Qj, Qj_d) );
            if(Qi_l != -1)
                noninteraction_edges.push_back( make_pair(Qi, Qi_l) );
            //4-1
            if( abs(yi-yj)==1  )
            {
                if(Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_r) );
                if(Qi_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_r) );
                if(Qj_r != -1 && Qj_2r != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_2r) );
                if(Qj_r != -1 && Qj_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qj_r, Qj_dr) );
                if(Qi_l != -1 && Qj_l != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qj_l) );
                if(Qi_r != -1 && Qj_r != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qj_r) );
                if(Qi_u != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u,   Qi_ul) );
                if(Qi_l != -1 && Qi_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qi_l, Qi_ul) );
            }
            //4-2
            else
            {
                if(Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qj, Qj_u) );
                if(Qi_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi, Qi_d) );
                if(Qj_u != -1 && Qj_2u != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_2u) );
                if(Qj_u != -1 && Qj_ul != -1)
                    noninteraction_edges.push_back( make_pair(Qj_u, Qj_ul) );
                if(Qi_u != -1 && Qj_u != -1)
                    noninteraction_edges.push_back( make_pair(Qi_u,   Qj_u) );
                if(Qi_d != -1 && Qj_d != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d, Qj_d) );
                if(Qi_r != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_r, Qi_dr) );
                if(Qi_d != -1 && Qi_dr != -1)
                    noninteraction_edges.push_back( make_pair(Qi_d,  Qi_dr) );
            } 
        }
    }

}

int Qcircuit::QMapper::return_layout_score(Graph& layout_graph, bool new_layout)
{

    int temp_layout_score = 0;
    int layout_score      = 0;
    vector<pair<int,int>> noninteraction_edges;
    for(auto& edge : layout_graph.edgeset)
    {
        temp_layout_score = 0;
        int e_id     = edge.second.getid();
        int e_source = edge.second.getsourceid(); //q_source
        int e_target = edge.second.gettargetid(); //q_target
        int Q_source, Q_target;
        if(!new_layout)
        {
            Q_source = layout_L[e_source];        //Q_source
            Q_target = layout_L[e_target];        //Q_target
        }
        else
        {
            Q_source = new_layout_L[e_source];        //Q_source
            Q_target = new_layout_L[e_target];        //Q_target
        }
        int weight   = edge.second.getweight();   //interaction weight
        temp_layout_score += weight;
        
        noninteraction_edges.clear();
        return_noninteraction_edges(noninteraction_edges, Q_source, Q_target);

        for(auto& pair : noninteraction_edges)
        {
            for(auto& e : layout_graph.edgeset)
            {
                if(!new_layout)
                {
                    if(pair.first == layout_L[e.second.getsourceid()] && pair.second == layout_L[e.second.gettargetid()])
                        temp_layout_score -= e.second.getweight();
                    else if(pair.second == layout_L[e.second.getsourceid()] && pair.first == layout_L[e.second.gettargetid()])
                        temp_layout_score -= e.second.getweight();
                }
                else
                {
                    if(pair.first == new_layout_L[e.second.getsourceid()] && pair.second == new_layout_L[e.second.gettargetid()])
                        temp_layout_score -= e.second.getweight();
                    else if(pair.second == new_layout_L[e.second.getsourceid()] && pair.first == new_layout_L[e.second.gettargetid()])
                        temp_layout_score -= e.second.getweight();
                }
            }
        }
        layout_score += temp_layout_score;
    }
    return layout_score;
}

bool Qcircuit::QMapper::make_layout(int move_x, int move_y, bool isReflected)
{
    //initiallize layout -> -1
    for(int i=0; i<positions; i++)
        new_layout_L[i] = -1; //logical -> physical
    for(int i=0; i<positions; i++)
        new_qubit_Q[i ] = -1; //physical -> logical
    
    //for(int i=0; i<nqubits; i++)
    for(int i=0; i<positions; i++)
        new_layout_L[i] = layout_L[i]; //logical -> physical
    for(int i=0; i<positions; i++)
        new_qubit_Q[i ] = qubit_Q[i ]; //physical -> logical

    //sliding for move_x & move_y
    for(int i=0; i<nqubits; i++)
    {
        int Qi = new_layout_L[i];
        int xi = p_to_xy[Qi].first;
        int yi = p_to_xy[Qi].second;
       
        //isReflected check
        if(isReflected)
        {
            xi *= -1;
            int new_Qi = xy_to_p(xi, yi);
            if(new_Qi == -1)
                return false;
        }
        
        int new_xi = xi + move_x;
        int new_yi = yi + move_y;
        int new_Qi = xy_to_p(new_xi, new_yi);
        if(new_Qi == -1)
            return false;
        
        //layout update
        new_layout_L[i    ] = new_Qi;
    }
    for(int i=0; i<positions; i++)
        new_qubit_Q[i] = -1;
    for(int i=0; i<positions; i++)
    {
        if(new_layout_L[i] != -1)
            new_qubit_Q[ new_layout_L[i] ] = i;
    }

    return true;

    
}
