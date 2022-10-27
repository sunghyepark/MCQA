//mappingPreparation.cpp

#include "circuit.h"

using namespace std;
using namespace Qcircuit;

bool compare_cost(pair<int, int> a, pair<int, int> b)
{
    return a.first < b.first;
}

bool compare_cost_d(pair<int, int> a, pair<int, int> b)
{
    return a.first > b.first;
}

bool compare_graph_bfs(pair< pair<int, int>, int> a, pair< pair<int, int>, int> b)
{
    if(a.first.first == b.first.first)
        return a.first.second < b.first.second;
    else
        return a.first.first < b.first.first;
}

bool compare_graph_bfs_number(pair< pair<int, int>, int> a, pair< pair<int, int>, int> b)
{
    if(a.first.first == b.first.first)
        return a.first.second > b.first.second;
    else
        return a.first.first < b.first.first;
}


bool compare(const pair<int, int>& a, const pair<int, int>& b)
{
    //if the first number is same
    if(a.first == b.first)
        return a.second < b.second;
    return a.first > b.first;
}

///////////////////////////////
void Qcircuit::QMapper::generate_bfs_queue(queue<int>& queue, Graph& graph, int start, bool order)
{
    vector< pair< pair<int, int>, int> >  bfs_sort_vector; // <dist, weight> , nodeid(table idx)
    for(int i=0; i<nqubits; i++)
    {
        if(layout_L[i]!=-1 && i!=start)
            continue;
        bfs_sort_vector.push_back(make_pair(make_pair(graph.dist[start][i], graph.bfs_weight(start, i)), i));
    }

    if(order) //order: small -> large
        sort(bfs_sort_vector.begin(), bfs_sort_vector.end(), compare_graph_bfs);
    else //number: large -> small
        sort(bfs_sort_vector.begin(), bfs_sort_vector.end(), compare_graph_bfs_number);
   
    queue.push(start);  //for graph matching 

    for(int i=0; i<bfs_sort_vector.size(); i++)
        queue.push(bfs_sort_vector[i].second);
    
    if(queue.front() == start) queue.pop();
}

void Qcircuit::QMapper::sort_degree(vector<int>& candi_loc, Graph& graph)
{
    vector<pair<int, int>> degree_candi_loc;
    vector<int> temp_candi_loc;
    int degree;
    int min_dist = 1;
    
    //** initial candi loc
    for(auto v : candi_loc)
    {

        degree = 0;
        for(int i = 0; i < graph.node_size; i++)
        {
            if(i == v) continue;
            if(graph.dist[v][i] > min_dist) continue;
            if(qubit_Q.count(i) && qubit_Q[i] == -1)
                degree++;
        }
        //cout << v << "(degree: " << degree << "), ";
        degree_candi_loc.push_back(make_pair(degree, v));
    }
    sort(degree_candi_loc.begin(), degree_candi_loc.end(), compare);
    for(auto v : degree_candi_loc)
    {
        temp_candi_loc.push_back(v.second);
    }
    candi_loc = temp_candi_loc;
}

int Qcircuit::QMapper::sort_degree_return(vector<int>& candi_loc, Graph& graph)
{
    vector<pair<int, int>> degree_candi_loc;
    vector<int> temp_candi_loc;
    int degree;
    int min_dist = 1;
    int return_degree = 0;

    //** initial candi loc
    for(auto v : candi_loc)
    {

        degree = 0;
        for(int i = 0; i < graph.node_size; i++)
        {
            if(i == v) continue;
            if(graph.dist[v][i] > min_dist) continue;
            if(qubit_Q.count(i) && qubit_Q[i] == -1)
                degree++;
        }
        //cout << v << "(degree: " << degree << "), ";
        degree_candi_loc.push_back(make_pair(degree, v));
        return_degree += degree;
    }
    
    sort(degree_candi_loc.begin(), degree_candi_loc.end(), compare);
    for(auto v : degree_candi_loc)
    {
        temp_candi_loc.push_back(v.second);
    }
    candi_loc = temp_candi_loc;

    return return_degree;
}



///////////////////////////////
void Qcircuit::QMapper::make_Glist(Circuit& dgraph)
{
    Glist.clear();
    for(int j = 0; j<positions; j++)
    {
        list<pair<int, int>> temp;
        Glist.push_back(temp);
    }
    for(const auto& gate : dgraph.nodeset)
    {
        if(gate.control != -1) // two-qubit gate
        {
            const int id      = gate.id;
            const int target  = gate.target;
            const int control = gate.control;
            Glist[target].push_back(make_pair(id, 2));
            Glist[control].push_back(make_pair(id, 2));
            Glist[target].push_back(make_pair(-1, 1));
            Glist[control].push_back(make_pair(-1, 1));
        }
        else // single qubit gate
        {
            const int id     = gate.id;
            const int target = gate.target; //control=-1
            Glist[target].push_back(make_pair(id, 1));
        }
    }
}

void Qcircuit::QMapper::make_CNOT(bool i)
{
    Dgraph_cnot.nodeset.clear();
    if(i)
    {
        //cout << "** CNOT set **" << endl;
        for(const auto& gate : Dgraph.nodeset)
        {
            if(gate.control == -1) continue;
            const int target  = gate.target;
            const int control = gate.control;
    
#if GATESET_IIC_JKU || GATESET_DAC2022
            add_cnot(control, target, Dgraph_cnot);
            add_cnot_num--;
#elif GATESET_QUTECH
            add_cz(control, target, Dgraph_cnot);
            add_cz_num--;
#endif
        }
    }
}

Graph Qcircuit::QMapper::generate_interaction_graph(bool i, int n)
{
    bool** gen_edge = new bool*[nqubits]; //edge generated (for interaction graph)
    for(int i = 0; i < nqubits; i++)
        gen_edge[i] = new bool[nqubits];

    //gen_edge initialize
    for(int i = 0; i < nqubits; i++)
        for(int j = 0; j < nqubits; j++)
            if(i == j)
                gen_edge[i][j] = true;
            else
                gen_edge[i][j] = false;

    Graph g;
    for(int i = 0; i < nqubits; i++)
        g.addnode();

    int dgraphSize = Dgraph_cnot.nodeset.size();
    int m = dgraphSize / n;
    int remainder = dgraphSize % n;

    vector<double> costWeightTable;
    double costParam = 0.9;
    double cost = 1.0;
    costWeightTable.clear();
    for(int i=0; i<n+1; i++)
    {
        costWeightTable.push_back(cost);
        cost *= costParam;
    }
    
    //cout << "** Generate MixGraph ** " << endl;
    int for_i=0;
    int for_n=0;
    double weight_cost;
    for(const auto& gate : Dgraph_cnot.nodeset)
    {
        if(for_i==m)
        {
            for_i=0;
            for_n++;
        }
        weight_cost = costWeightTable[for_n];
        const int id = gate.id;
        const int target = gate.target;
        const int control = gate.control;
        if(gen_edge[target][control])
        {
            for(auto& e : g.edgeset)
            {
                if(e.second.gettargetid()==target && e.second.getsourceid()==control)
                {
                    int weight = e.second.getweight();
                    e.second.setweight(weight+weight_cost*10);
                }
                else if(e.second.gettargetid()==control && e.second.getsourceid()==target)
                {
                    int weight = e.second.getweight();
                    e.second.setweight(weight+weight_cost*10);
                }
                else continue;
            }
        }
        else
        g.addedge(target, control, weight_cost*10);
        gen_edge[target][control] = true;
        gen_edge[control][target] = true;

        for_i++;
    }

    for(int i = 0; i < nqubits; i++)
        delete[] gen_edge[i];
    delete[] gen_edge;

    return g;

}


Graph Qcircuit::QMapper::make_layoutGraph(bool i, bool new_layout)
{
    bool** gen_edge = new bool*[nqubits]; //edge generated
    for(int i=0; i<nqubits; i++)
        gen_edge[i] = new bool[nqubits ];

    //gen_edge initialize
    for(int i=0; i<nqubits; i++)
        for(int j=0; j<nqubits; j++)
            if(i==j)
                gen_edge[i][j] = true ;
            else
                gen_edge[i][j] = false;

    Graph g;
    for(int i=0; i<nqubits; i++)
        g.addnode();

    for(int i=0; i<nqubits; i++)
    {
        int qi, pi;
        qi = i;
        if(new_layout)
            pi = new_layout_L[qi];
        else
            pi = layout_L[qi];

        for(int j=0; j<nqubits; j++)
        {
            if(i!=j && !gen_edge[i][j])
            {
                int qj, pj;
                int weight;
                qj = j;
                if(new_layout)
                    pj = new_layout_L[qj];
                else
                    pj = layout_L[qj];
                if(coupling_graph.dist[pi][pj] == 1)
                {
                    weight = get_interaction_weight(qi, qj);
                    g.addedge(qi, qj, weight);
                    gen_edge[qi][qj] = true;
                    gen_edge[qj][qi] = true;
                }
            }
        }
    }

    for(int i=0; i<nqubits; i++)
        delete[] gen_edge[i];
    delete[] gen_edge;

    return g;
}


void Qcircuit::QMapper::make_ref_loc(vector<int>& ref_loc, Graph& graph, int start, bool order)
{
    vector<pair<int, int>> temp; //edge weight, id
    //convert table idx to graph node id
    int g_id = graph.nodeid[start];
    Node& node = graph.nodeset[g_id];
    for(auto e_id : node.edges)
    {
        Edge& edge = graph.edgeset[e_id];
        int id = (edge.source->getid() == start) ? edge.target->getid() : edge.source->getid();
        if(layout_L[id] != -1) temp.push_back(make_pair(edge.getweight(), layout_L[id]));
    }
    if(order)
        sort(temp.begin(), temp.end(), compare_cost);
    else
        sort(temp.begin(), temp.end(), compare_cost_d);

    for(auto v : temp)
        ref_loc.push_back(v.second);
    /*
    cout << " ref_loc: ";
    for(auto v : temp)
    { 
        cout << v.second << ", ";
    }
    cout << endl;
    */
}

void Qcircuit::QMapper::make_candi_loc(int current_qc, vector<int>& candi_loc, vector<int>& ref_loc, Graph& coupling_graph, Graph& interaction_graph, bool equal_order)
{
    vector<pair<int, int>> distance_candi_loc;
    vector<int> initial_candi_loc;

    //initial candi_loc
    for(int i = 0; i < coupling_graph.node_size; i++)
    {
        if(qubit_Q.count(i) && qubit_Q[i] == -1)
            initial_candi_loc.push_back(i);
    }
    //cout << "initial candi loc: " << initial_candi_loc << endl;
        //sub_num --> candi_loc
    if(equal_order)
    {
        for(auto r : ref_loc)
        {
            for(auto v : initial_candi_loc)
            {
                int dist = coupling_graph.dist[v][r];
                distance_candi_loc.push_back(make_pair(dist, v));
            }
            sort(distance_candi_loc.begin(), distance_candi_loc.end(), compare_cost);
            initial_candi_loc.clear();
            for(auto v : distance_candi_loc)
            {
                if(v.first == distance_candi_loc[0].first)
                    initial_candi_loc.push_back(v.second);
            }
            distance_candi_loc.clear();
            //cout << "initial candi loc: " << initial_candi_loc << endl;
        }
    }
    //dist_sum --> candi_loc
    else
    {
        for(auto v : initial_candi_loc)
        {
            int dist_sum = 0;
            int dist;
            for(auto r : ref_loc)
            {
                dist = coupling_graph.dist[v][r];
                dist_sum += dist;
            }
            distance_candi_loc.push_back(make_pair(dist_sum, v));
        }
        sort(distance_candi_loc.begin(), distance_candi_loc.end(), compare_cost);
        initial_candi_loc.clear();
        for(auto v : distance_candi_loc)
        {
            if(v.first == distance_candi_loc[0].first)
                initial_candi_loc.push_back(v.second);
        }
        distance_candi_loc.clear();
        //cout << "initial candi loc: " << initial_candi_loc << endl;
    }
    candi_loc = initial_candi_loc;
}

