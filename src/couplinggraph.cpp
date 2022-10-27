//couplinggraph.cpp

#include "circuit.h"
using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::select_coupling_graph(ARCHITECTURE archi)
{
    switch(archi)
    {
        case ARCHITECTURE::Surface_17:
            build_graph_Surface_17();
            break;
        default:
            cout << "No architecture specified!" << endl;
            exit(1);
    }
}

void Qcircuit::QMapper::build_graph_Surface_17()
{    
    positions = 17;
    //add node
    coupling_graph.addnode(2); //0
    coupling_graph.addnode(1); //1
    coupling_graph.addnode(1); //2
    coupling_graph.addnode(1); //3
    coupling_graph.addnode(2); //4
    coupling_graph.addnode(2); //5
    coupling_graph.addnode(2); //6
    coupling_graph.addnode(3); //7
    coupling_graph.addnode(3); //8
    coupling_graph.addnode(3); //9
    coupling_graph.addnode(2); //10
    coupling_graph.addnode(2); //11
    coupling_graph.addnode(2); //12
    coupling_graph.addnode(1); //13
    coupling_graph.addnode(1); //14
    coupling_graph.addnode(1); //15
    coupling_graph.addnode(2); //16

    //addedge
    coupling_graph.addedge(0, 2);
    coupling_graph.addedge(0, 3);
    coupling_graph.addedge(1, 4);
    coupling_graph.addedge(1, 5);
    coupling_graph.addedge(2, 5);
    coupling_graph.addedge(2, 6);
    coupling_graph.addedge(3, 6);
    coupling_graph.addedge(4, 7);
    coupling_graph.addedge(5, 7);
    coupling_graph.addedge(5, 8);
    coupling_graph.addedge(6, 8);
    coupling_graph.addedge(6, 9);
    coupling_graph.addedge(7, 10);
    coupling_graph.addedge(8, 10);
    coupling_graph.addedge(8, 11);
    coupling_graph.addedge(9, 11);
    coupling_graph.addedge(9, 12);
    coupling_graph.addedge(10, 13);
    coupling_graph.addedge(10, 14);
    coupling_graph.addedge(11, 14);
    coupling_graph.addedge(11, 15);
    coupling_graph.addedge(12, 15);
    coupling_graph.addedge(13, 16);
    coupling_graph.addedge(14, 16);
    cout <<  " # This graph is: Surface 17" << endl;

    //make same_frequency group
    vector<int> f0;
    vector<int> f1 = {1,2,3,13,14,15};
    vector<int> f2 = {0,4,5,6,10,11,12,16};
    vector<int> f3 = {7,8,9};
    same_freq_group.push_back(f0);
    same_freq_group.push_back(f1);
    same_freq_group.push_back(f2);
    same_freq_group.push_back(f3);
}

void Qcircuit::QMapper::build_coordinate_Surface_17()
{
    p_to_xy.clear();
    p_to_xy.push_back( make_pair(  2,  1) ); //Q0
    p_to_xy.push_back( make_pair(  0,  2) ); //Q1
    p_to_xy.push_back( make_pair(  1,  1) ); //Q2
    p_to_xy.push_back( make_pair(  2,  0) ); //Q3
    p_to_xy.push_back( make_pair( -1,  2) ); //Q4
    p_to_xy.push_back( make_pair(  0,  1) ); //Q5
    p_to_xy.push_back( make_pair(  1,  0) ); //Q6
    p_to_xy.push_back( make_pair( -1,  1) ); //Q7
    p_to_xy.push_back( make_pair(  0,  0) ); //Q8
    p_to_xy.push_back( make_pair(  1, -1) ); //Q9
    p_to_xy.push_back( make_pair( -1,  0) ); //Q10
    p_to_xy.push_back( make_pair(  0, -1) ); //Q11
    p_to_xy.push_back( make_pair(  1, -2) ); //Q12
    p_to_xy.push_back( make_pair( -2,  0) ); //Q13
    p_to_xy.push_back( make_pair( -1, -1) ); //Q14
    p_to_xy.push_back( make_pair(  0, -2) ); //Q15
    p_to_xy.push_back( make_pair( -2, -1) ); //Q16
}

int Qcircuit::QMapper::xy_to_p(int x, int y)
{
    if(x==-2)
    {
        if(y==0)
            return 13;
        else if(y==-1)
            return 16;
        else return -1;
    }
    else if(x==-1)
    {
        if(y==2)
            return 4;
        else if(y==1)
            return 7;
        else if(y==0)
            return 10;
        else if(y==-1)
            return 14;
        else return -1;
    }
    else if(x==0)
    {
        if(y==2)
            return 1;
        else if(y==1)
            return 5;
        else if(y==0)
            return 8;
        else if(y==-1)
            return 11;
        else if(y==-2)
            return 15;
        else return -1;
    }
    else if(x==1)
    {
        if(y==1)
            return 2;
        else if(y==0)
            return 6;
        else if(y==-1)
            return 9;
        else if(y==-2)
            return 12;
        else return -1;

    }
    else if(x==2)
    {
        if(y==1)
            return 0;
        else if(y==0)
            return 3;
        else return -1;
    }
    else return -1;
}
