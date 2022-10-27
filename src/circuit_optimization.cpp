//circuit_optimization.cpp
//added by dhkim, 2022.04.28

#include "circuit.h"


namespace Qcircuit {

using namespace std;

void QMapper::circuit_optimization(Circuit& dgraph, bool is_initial) {
    
    enum class Axis {
        NONE, X, Y
    };

    class RotationInfo {
    public:
        double get_angle() { return this->angle; }
        bool is_axis(Axis axis) { return this->axis == axis; }

        void set_angle(double angle) { this->angle = normalize_angle(angle); }
        void set_axis(Axis axis) { this->axis = axis; }

        void add_angle(double angle) { set_angle(this->angle + angle); }
        
    private:
        double normalize_angle(double x) {
            // noramalize angle to [-pi, pi) 
            return atan2(sin(x), cos(x));
        }

        double angle{ 0 };
        Axis axis{ Axis::NONE };
    };

    auto gatetype_to_axis = [](GATETYPE g) {
        switch (g) {
            case GATETYPE::RX: return Axis::X;
            case GATETYPE::RY: return Axis::Y;
            default: return Axis::NONE;
        }
    };

    auto create_gate = [](RotationInfo &info, int q) {
        Gate g;
        g.target = q;
        g.control = -1;
        g.angle = info.get_angle();
        
        if (info.is_axis(Axis::X)) {
            g.type = GATETYPE::RX;
        } else if (info.is_axis(Axis::Y)) {
            g.type = GATETYPE::RY;
        } else {
            cout << "warning: create_gate() - something went wrong..." << endl;
        }
        snprintf ( g.output_type, 127, "(%.2f)", g.angle);
        return g;
    };

    // create new nodeset
    vector<Gate> new_nodeset;
    vector<RotationInfo> rotation(is_initial ? nqubits : positions);

    for (const Gate &gate : dgraph.nodeset) {
        RotationInfo &info = rotation.at(gate.target);
        Axis axis = gatetype_to_axis(gate.type);
        if (info.is_axis(Axis::NONE)) {
            if (axis != Axis::NONE) { 
                // axis change: NONE->X / NONE->Y
                info.set_axis(axis);
                info.set_angle(gate.angle);
            }
        } else if (info.is_axis(axis)) {
            // no axis change: X->X / Y->Y
            info.add_angle(gate.angle);
        } else if (axis == Axis::NONE) {
            // axis change: X->NONE / Y->NONE
            RotationInfo &info_c = rotation.at(gate.control);
            if (info_c.get_angle() != 0) {
                new_nodeset.push_back(create_gate(info_c, gate.control));
            }
            if (info.get_angle() != 0) {
                new_nodeset.push_back(create_gate(info, gate.target));
            }         
            info_c.set_axis(axis);
            info_c.set_angle(0);
            info.set_axis(axis);
            info.set_angle(0);
            new_nodeset.push_back(gate);    
        } else {
            // axis change: X->Y / Y->X
            if (info.get_angle() != 0) {
                new_nodeset.push_back(create_gate(info, gate.target));
            }
            info.set_axis(axis);
            info.set_angle(gate.angle);
        } 
    }

    // add leftover single-qubit gates to new_nodeset
    for (int q = 0; q < nqubits; q++) {
        RotationInfo &info = rotation[q];
        if (info.is_axis(Axis::NONE)) 
            continue;

        if (info.get_angle() != 0) {
            new_nodeset.push_back(create_gate(info, q));
        }
    }

    // re-index the gate id
    for (int i = 0; i < new_nodeset.size(); i++) {
        new_nodeset[i].id = i;
    }
    
    // replace original nodeset with new_nodeset
    dgraph.nodeset = new_nodeset;
}

} // namespace Qcircuit
