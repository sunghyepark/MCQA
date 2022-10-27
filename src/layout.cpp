//layout.cpp
//added by dhkim, 2022.04.22

#include "circuit.h"

namespace Qcircuit {

Layout::Layout() { }

Layout::Layout(std::map<int, int> &layout_map) {
    n_physical = layout_map.size();
    log2phy_vec = std::vector<int>(n_physical, -1);
    phy2log_vec = std::vector<int>(n_physical, -1);

    for (const auto& q : layout_map) {
        log2phy_vec.at(q.first) = q.second;
        phy2log_vec.at(q.second) = q.first;
    }

}

void Layout::swap_physical(const int p1, const int p2) {
    const int q1 = log2phy_vec.at(p1);
    const int q2 = log2phy_vec.at(p2);

    std::swap(log2phy_vec.at(q1), log2phy_vec.at(q2));
    std::swap(phy2log_vec.at(p1), phy2log_vec.at(p2));
}

void Layout::swap_logical(const int q1, const int q2) {
    const int p1 = log2phy_vec.at(q1);
    const int p2 = log2phy_vec.at(q2);

    std::swap(log2phy_vec.at(q1), log2phy_vec.at(q2));
    std::swap(phy2log_vec.at(p1), phy2log_vec.at(p2));
}

int Layout::phy2log(const int physical) {
    if (0 <= physical && physical < n_physical)
        return phy2log_vec.at(physical);
    return -1;
}

int Layout::log2phy(const int logical) {
    if (0 <= logical && logical < n_physical)
        return log2phy_vec.at(logical);
    return -1;
}

void Layout::print_log2phy() {
    for (int i = 0; i < n_physical; i++)
        std::cout << " q" << setw(2) << i << "     ->     p" << setw(2) << log2phy(i) << std::endl;
    std::cout << std::endl;
}

void Layout::print_phy2log() {
    for (int i = 0; i < n_physical; i++)
        std::cout << " p" << setw(2) << i << "     ->     q" << setw(2) << phy2log(i) << std::endl;
    std::cout << std::endl;
}

} // namespace Qcircuit
