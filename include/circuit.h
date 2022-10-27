//circuit.h

/* --------------------------------------- */
/* This script is written by Sunghye Park  */
/* 2022.03.21                              */
/*                                         */
/* shpark96@postech.ac.kr                  */
/* --------------------------------------- */
#ifndef __CIRCUIT__
#define __CIRCUIT__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <set>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <iterator>
#include <utility>

// NOTE: added gateset option (dhkim, 2022.04.16)
#define GATESET_IIC_JKU 0
#define GATESET_DAC2022 0
#define GATESET_QUTECH  1

// NOTE: added macro for PI (dhkim, 2022.04.19)
#define QASM_PI 3.141592653589793238463

using namespace std;

namespace Qcircuit
{
    //initial mapping
    enum class INITIALTYPE
    {
        IDENTITY_MAPPING,
        RANDOM_MAPPING,
        GRAPH_MATCHING_PROPOSED,
    };

    //coupling graph
    enum class ARCHITECTURE
    {
        Surface_17,
    };

    enum class GATETYPE
    {
#if GATESET_IIC_JKU 
// >>>>>>>>>>>>>>>>>>>> IIC-JKU Original gateset >>>>>>>>>>>>>>>>>>>>
        U, CNOT, 
// <<<<<<<<<<<<<<<<<<<< IIC-JKU Original gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_DAC2022
// >>>>>>>>>>>>>>>>>>>> DAC2022 target gateset >>>>>>>>>>>>>>>>>>>>
        X, Y, Z, H, S, SDG, T, TDG, RX, RZ, CNOT
// <<<<<<<<<<<<<<<<<<<< DAC2022 target gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_QUTECH
// >>>>>>>>>>>>>>>>>>>> QuTech gateset >>>>>>>>>>>>>>>>>>>>
        RX, RY, CZ
// <<<<<<<<<<<<<<<<<<<< QuTech gateset <<<<<<<<<<<<<<<<<<<<
#endif
    };


    //Circuit info (node, edge)
    struct Gate
    {
        int id; 
        int control, target;                //logical qubits //if single qubit control = -1
        GATETYPE type;
            
        //// use just in QASM // double theta, double phi, double lambda //
        char output_type[128];

        // NOTE: added by dhkim, 2022.04.28
        double angle{ 0 }; // for circuit_optimization.cpp
    };
  
    struct Circuit
    {
        vector<Gate> nodeset;   /* gate set */
    };

    struct Node;
    struct Edge;

    struct Node
    {
        private:
        ////variables
        static int global_id;
        int id;
        int weight;
        public:
        set<int> edges;

        public:
        ////constructors
        Node();
        Node(int weight);

        ////getters
        int getglobalid();
        int getid();
        int getweight();
        void setweight(int weight);

        static void init_global_id();
    };
    
    struct Edge
    {
        private:
        ////variables
        static int global_id;
        int id;
        int weight;
        public:
        Node* source;
        Node* target;

        public:
        ////constructors
        Edge();
        Edge(Node* source, Node* target, int weight);

        ////getters
        int getglobalid();
        int getid();
        int getweight();
        int getsourceid();
        int gettargetid();
        void setweight(int weight);
        void setsource(Node* source);
        void settarget(Node* target);

        static void init_global_id();
    };

    struct Graph
    {
        ////graph node and edge
        map<int, Node> nodeset;
        map<int, Edge> edgeset;
        
        map<int, int> nodeid;       //table index <--> graph node id
        int node_size;
        int** dist;

        int center;

        public:
        ////constructors
        Graph();
        
        //methods
        //graph.cpp
        void init_global_id();
        int addnode(int weight = 1);
        bool deletenode(int id);
        bool addedge(int source, int target, int weight = 1);
        bool deleteedge(int id);

        //graphFunction.cpp
        void build_dist_table();
        void print_dist_table();
        void delete_dist_table();
        
        int bfs(int start, int end);
        int bfs_weight(int start, int end);
        
        int generate_graph_center();
        pair<int, bool> graph_characteristics(bool i);
    };

    class Layout {
        // NOTE: added by dhkim, 2022.04.22
    public:
        Layout();
        Layout(std::map<int, int> &layout_map);
        
        void swap_logical(const int q1, const int q2);
        void swap_physical(const int p1, const int p2);

        int phy2log(const int q_logical);
        int log2phy(const int q_physical);

        void print_log2phy();
        void print_phy2log();

    private:
        int n_physical;
        std::vector<int> log2phy_vec;
        std::vector<int> phy2log_vec;
    };

    class QMapper
    {
        public:
        //////////// variables
            //file names
            string fileName_input;
            string fileName_output;
            
            //cost parameters
            double param_alpha; //discount factor when calculating swap cost
            double param_beta; //select between BRIDGE & (SWAP or MOV)
            vector<double> cost_weight_table;

            //quantity of logical, physical qubits
            unsigned int nqubits;
            unsigned int positions;
            
            //Quantum circuit
            Circuit Dgraph;
            Circuit Dgraph_cnot;

                //-> for initial mapping
            Graph interactionGraph; //quantum circuit -> interaction graph
            Graph layoutGraph;      //graph matching (interactionGraph -> coupling_graph)
            Graph new_layoutGraph;  //frequency matching

            //Coupling graph
            Graph coupling_graph;
            vector< vector<int> > same_freq_group; //frequency_coupling_graph[f_i] = {Q_a, ... Q_z} (same frequency group)
            
            //Quantum circuit -> Data structure (for main mapping)
            vector<list<pair<int, int>> > Glist;
            
            //initial mapping
            map<int, int> layout_L; //L[logical] = physical
            map<int, int> qubit_Q;  //Q[physical] = logical
            map<int, int> new_layout_L; //for freq. matching
            map<int, int> new_qubit_Q;  //for freq. matching
            Layout layout;
                //for final latency calculation
            Layout initial_layout; 
                //for freq matching
            vector< pair<int, int> > p_to_xy; //Qi -> (x, y)
            vector<int> freq_mod; //[f3, f2, f1, f2]

            //Final quantum circuit
            Circuit FinalCircuit; //final quantum circuit

            //from mapping.cpp
            int node_id;
           
#if GATESET_IIC_JKU || GATESET_DAC2022
            int add_cnot_num;
#elif GATESET_QUTECH
            int add_cz_num;
#endif
            int add_swap_num;
            int add_mov_num;
            int add_bridge_num;
            int latency;

        /////////// functions
            //parser.cpp
            void parsing(int argc, char** argv);
            void QASM_parse(string filename);
            
            //couplinggraph.cpp
            void select_coupling_graph(ARCHITECTURE archi);
            void build_graph_Surface_17();
            void build_coordinate_Surface_17();
            int xy_to_p(int x, int y);
            
            //circuit_optimization.cpp
            // NOTE: added by dhkim, 2022.04.28
            void circuit_optimization(Circuit& dgraph, bool is_initial);

            //initial_mapping.cpp
            void initial_mapping(ARCHITECTURE archi, bool prepro, bool BRIDGE_MODE, INITIALTYPE initial_type);
            void graph_matching_processing_proposed(ARCHITECTURE archi);
            void identical_mapping();
            void random_mapping();

            int get_interaction_weight(int sourceid, int targetid);
            int generate_candi_layout(vector<int>& candi_loc1, vector<int>& candi_loc2, bool dist_max, bool forward);
            void print_layout(const map<int, int> &layout);
            void print_qubit(const map<int, int> &qubit);
            void print_layout();
            void print_qubit();

            int return_freq_idx(int x, int y);
            void return_noninteraction_edges(vector<pair<int,int>>& noninteraction_edges, int pi, int pj); //pi(xi,yi), p2(xj,yj) -> return noninteraction_edges
            int return_layout_score(Graph& layout_graph, bool new_layout);
            bool make_layout(int move_x, int move_y, bool isReflected);

            //mapping.cpp
            void Circuit_Mapping(Circuit& dgraph, ARCHITECTURE archi,  bool BRIDGE_MODE);
            void directly_execute_gate(list<int>& g_scheduled, list<int>& g_sat, vector<int>& frozen, vector<int>& frozen_duration, vector<bool>& updated_Glist, Circuit& dgraph);
                //update gate list
            void update_Glist(vector<int>& frozen_duration, vector<bool>& updated_Glist, Circuit& dgraph);
            void update_front_n_pair_gates(list<int>& g_head, list<int>& g_pair, list<int>& g_sat, list<int>& g_vio, vector<int>& frozen);
            bool check_connectivity(list<int>& g_head, list<int>& g_pair, list<int>& g_sat, list<int>& g_vio, vector<int>& frozen, Circuit& dgraph);
            void update_bridge_gates(list<int>& g_bridge, list<int>& g_vio, vector<int>& frozen_constraints, Circuit& dgraph);
                //constraint resolve
            void generate_swap_candi_list(list<int>& g_vio, vector< pair< pair<int, int>, int> >& swap_candi_list, vector<int>& frozen_constraints, Circuit& dgraph);
            void generate_mov_candi_list(list<int>& g_vio, vector< pair< pair<int, int>, int> >& mov_candi_list, vector<int>& frozen_constraints, Circuit& dgraph);
            bool check_connected_graph(int p_original, int p_new);
                //gate scheduling
            void update_frozen(vector<int>& frozen);
            void sort_priority_gates(vector< pair<pair<int, int>, pair<int, int>> >& g_temp, list<int>& g_scheduled, vector<int>& frozen_duration);
            void generate_scheduled_duration(list<int>& gates_list, vector< pair<pair<int, int>, pair<int, int>> >& g_temp, vector<int>& frozen_duration, Circuit& dgraph);
            void update_frozen_duration(vector<int>& frozen_duration);
            void generate_freq_same(vector<vector<int>>& same_freq_qubits); //[m]: physical qubits -> same_freq. groups
            void generate_freq_relations(vector< pair<vector<int>, vector<int>> >& freq_relations); //[m]: physical qubits -> first: lower freq. groups / second: higher freq. groups
            void generate_scheduled_frequency(list<int>& gates_list, vector<int>& frozen_frequency, Circuit& dgraph);
            void update_frozen_frequency(vector<int>& frozen_frequency);
                //main mapping code organize
            void print_gates(list<int>& gates, int gate_type);              //gate_type) 0: front, 1: pair, 2: satisfied, 3: violated, 4: bridge, 5: scheduled
            void print_temp_scheduled_gates(vector<pair<pair<int,int>, pair<int,int>> >& g_temp);
            void print_gates_spec(list<int>& gates, Circuit& dgraph);
            void print_frozen(vector<int>& frozen_vector, int frozen_type); //frozen_type) 0: frozen, 1: frozen_duration, 2: frozen_frequency, 3: frozen_constraints;
            void print_Glist();
            void print_candi_list(vector< pair< pair<int,int>, int> >& candi_list, int candi_type); //candi_type) 0: swap, 1: mov
            void print_addGates(int gateid, int q_control, int q_target, int p_control, int p_target);
            void print_candi_pair(pair<int,int>& Pair, int gateid, int max_cost, int candi_type); //candi_type) 0: swap, 1: mov

            //mappingPreparation.cpp
            void make_Glist(Circuit& dgraph);
            void make_CNOT(bool i);
            Graph generate_interaction_graph(bool i, int n);
            Graph make_layoutGraph(bool i, bool new_layout); //new_layout=false -> using layout_L // true- -> using new_layout_L
            void sort_degree(vector<int>& candi_loc, Graph& graph);
            int sort_degree_return(vector<int>& candi_loc, Graph& graph);
            void generate_bfs_queue(queue<int> & queue, Graph& graph, int start, bool order);
            void make_ref_loc(vector<int>& ref_loc, Graph& graph, int start, bool order);
            void make_candi_loc(int current_qc, vector<int>& candi_loc, vector<int>& ref_loc, Graph& coupling_graph, Graph& interaction_graph, bool equal_order);

            //mappingFunction.cpp
            int cal_effect(const int q1, const int q2, const int p1, const int p2);
            double cal_MCPE(const pair<int, int> p, Circuit& dgraph, bool isswap); //isswap:1 -> swap/ isswap:0 -> mov
            void find_max_cost(pair<int,int>& SWAP, int& gateid, double& max_swap_cost,
                                                        vector< pair< pair<int,int>, pair<int,double> > >& v, list<int>& g_vio,
                                                        vector< pair< pair<int,int>, pair<int,double> > >& history);
            void add_cnot(int c_qubit, int t_qubit, Circuit& graph);
#if GATESET_QUTECH
            void add_cz(int c_qubit, int t_qubit, Circuit& graph);
#endif
            void add_swap(int q1, int q2, Circuit& graph);
            void add_mov(int q1, int q2, Circuit& graph);
            void add_bridge(int ps, int pt, vector<int>& frozen_constraints, Circuit& graph);
            void add_gate_to_Glist(int gate_mode, int update_nodeid, int q_control, int q_target, int q_middle);
            void _add_gate_to_Glist(int q1, int q2, int& gateid);

            //latency.cpp
            // NOTE: added by dhkim, 2022.05.03
            void update_latency(Circuit& dgraph, ARCHITECTURE archi);

            //outputwriter.cpp
            void FinalCircuit_info(Circuit& graph);
            void write_output(Circuit& graph);
    };
}

//----------------------------- for debug ----------------------------//
using namespace Qcircuit;

inline ostream& operator << (ostream& os, Node n)
{
    os << "Node[" << setw(3) << n.getid() << "] | weight: " << setw(2) << n.getweight();
    if(!n.edges.empty())
    {
        os << " | edge ids: ";
        for(auto e : n.edges)
            os << setw(3) << e << " ";
    }
    return os;
}

inline ostream& operator << (ostream& os, Edge e)
{
    os << "Edge[" << setw(3) << e.getid() << "] | weight: " << setw(2) << e.getweight();
    os << " | " << setw(3) << e.getsourceid() << " <--> " << setw(3) << e.gettargetid();
    return os;
}

inline ostream& operator << (ostream& os, map<int, Node> nodeset)
{
    for(auto& kv : nodeset)
        os << kv.second << endl;
    return os;
}

inline ostream& operator << (ostream& os, map<int, Edge> edgeset)
{
    for(auto& kv : edgeset)
        os << kv.second << endl;
    return os;
}

inline ostream& operator << (ostream& os, Graph graph)
{
    os << "||Graph||" << endl;
    os << "||Node Set||" << endl;
    os << graph.nodeset << endl;
    os << "||Edge Set||" << endl;
    os << graph.edgeset << endl;
    return os;
}

inline ostream& operator << (ostream& os, const GATETYPE& type)
{
    switch(type)
    {
#if GATESET_IIC_JKU
// >>>>>>>>>>>>>>>>>>>> IIC-JKU Original gateset >>>>>>>>>>>>>>>>>>>>
        case GATETYPE::U    : return os << "u"  ;
        case GATETYPE::CNOT : return os << "cx" ;
// <<<<<<<<<<<<<<<<<<<< IIC-JKU Original gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_DAC2022
// >>>>>>>>>>>>>>>>>>>> DAC2022 target gateset >>>>>>>>>>>>>>>>>>>>
        case GATETYPE::X    : return os << "x"  ;
        case GATETYPE::Y    : return os << "y"  ;
        case GATETYPE::Z    : return os << "z"  ;
        case GATETYPE::H    : return os << "h"  ;
        case GATETYPE::S    : return os << "s"  ;
        case GATETYPE::SDG  : return os << "sdg";
        case GATETYPE::T    : return os << "t"  ;
        case GATETYPE::TDG  : return os << "tdg";
        case GATETYPE::RX   : return os << "rx" ;
        case GATETYPE::RZ   : return os << "rz" ;
        case GATETYPE::CNOT : return os << "cx" ;
// <<<<<<<<<<<<<<<<<<<< DAC2022 target gateset <<<<<<<<<<<<<<<<<<<<
#elif GATESET_QUTECH
// >>>>>>>>>>>>>>>>>>>> QuTech gateset >>>>>>>>>>>>>>>>>>>>
        case GATETYPE::RX   : return os << "rx";
        case GATETYPE::RY   : return os << "ry";
        case GATETYPE::CZ   : return os << "cz";
// <<<<<<<<<<<<<<<<<<<< QuTech gateset <<<<<<<<<<<<<<<<<<<<
#endif // GATESET
        default             : return os;
    }
}

inline ostream& operator << (ostream& os, const Gate& gate)
{
    os << "Gate";
    if(gate.control == -1) //U, RX, RZ gate
        os << "(" << setw(5) << gate.id << ")\t" << setw(2) << gate.type << " " << "q[" << gate.target << "]";
    else //CNOT gate
        os << "(" << setw(5) << gate.id << ")\t" << setw(2) << gate.type << " " << "q[" << gate.control << "], q[" << gate.target << "]";
    return os;   
}

inline ostream& operator << (ostream& os, const vector<Gate>& nodeset)
{
    for(auto gate : nodeset)
        os << gate << endl;
    return os;
}

inline ostream& operator << (ostream& os, const queue<int> q)
{
    queue<int> test = q;
    while(!test.empty())
    {
        os << test.front() << " ";
        test.pop();
    }
    return os;
}

inline ostream& operator << (ostream& os, const list<int>& list)
{
    for(auto v : list)
        os << v << " ";
    return os;
}

inline ostream& operator << (ostream& os, const vector<list<int>>& list)
{
    for(int n=0; n<list.size(); n++)
    {
        os << "n[" << setw(4) << n << "]: " << "size: " << list[n].size() << " ";
        for(auto ll : list[n])
            os << "-> " << ll;
        os << endl;
    }
    return os;
}

inline ostream& operator << (ostream& os, const vector<list<pair<int, int>> >& list)
{
    for(int n=0; n<list.size(); n++)
    {
        os << "n[" << setw(4) << n << "]: " << "size: " << setw(4) << list[n].size() << " ";
        std::list<pair<int, int> >::const_iterator it = list[n].begin();
        int while_num = 15;
        while(while_num !=0 && it != list[n].end())
        {
            //os << "-> " << ll;
            if(it->first != -1)
                os << "-> (" << setw(3) << it->first << ", " << it->second << ")";
            else
                os << "-> - - - - ";
            it++;
            while_num--;
        }
        os << endl;
    }
    return os;
}

inline ostream& operator << (ostream& os, const vector<int>& list)
{
    for(auto v : list)
        os << v << " ";
    return os;
}

inline ostream& operator << (ostream& os, const vector<bool>& list)
{
    for(auto v : list)
        os << v << " ";
    return os;
}

inline ostream& operator << (ostream& os, const map<int, int>& m)
{
   for(auto kv : m)
        os << "map[" << setw(3) << kv.first  << "] = " << kv.second << endl;
    return os;
}

inline ostream& operator << (ostream& os, const vector<pair<int, int> >& llist)
{
    os << "{ ";
    for(auto v : llist)
        os << "(" << v.first << ", " << v.second << ") ";
    os << "}";
    return os;
}

inline ostream& operator << (ostream& os, const vector<pair<bool, bool> >& llist)
{
    os << "{ ";
    for(auto v : llist)
        os << "(" << v.first << ", " << v.second << ") ";
    os << "}";
    return os;
}

inline ostream& operator << (ostream& os, const vector<vector<int>>& list)
{
    for(auto v : list)
        os << v << " ";
    return os;
}

inline ostream& operator << (ostream& os, const pair<int, int> p)
{
    return os << "pair: " << p.first << " " << p.second;
}

inline ostream& operator << (ostream& os, const vector< pair< pair<int, int>, int> > v)
{
    for(auto kv : v)
        os << "SWAP: " << kv.first.first << " " << kv.first.second << " g_vio_id: " << kv.second << endl;
    return os;
}

inline ostream& operator << (ostream& os, const vector< pair< pair<int,int>, pair<int, int>> >  v)
{
        for(auto kv : v)
        {
            os << "gateid: " << kv.first.first << " logical_qubit:" << kv.first.second;
            os << " qubit_duration:" << kv.second.first << " Glist_length:" << kv.second.second << endl;
        }
        return os;
}

inline ostream& operator << (ostream& os, const vector< pair< pair<int,int>, pair<int, double>> >  v)
{
        for(auto kv : v)
        {
            os << "SWAP: " << kv.first.first << " " << kv.first.second;
            os << " g_vio_id: " << kv.second.first << " MCPE cost: " << kv.second.second << endl;
        }
        return os;
}
#endif // __CIRCUIT__
