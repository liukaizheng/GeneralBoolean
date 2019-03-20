#ifndef CELL_ADJACENCY_H
#define CELL_ADJACENCY_H

#include <tuple>
#include <Eigen/Dense>
#include <vector>
#include <set>

namespace bo {
template <typename DerivedC>
inline void cellAdjacency(
        const Eigen::MatrixBase<DerivedC> &per_patch_cells,
        const size_t & num_cells,
        std::vector<std::set<std::tuple<typename DerivedC::Scalar, bool, size_t>>> &adjacency_list)
{
    const size_t &num_patches = per_patch_cells.rows();
    adjacency_list.resize(num_cells);
    for(size_t i = 0; i < num_patches; i++) {
        const auto &pos = per_patch_cells(i, 0);
        const auto &neg = per_patch_cells(i, 1);
        adjacency_list[pos].emplace(neg, false, i);
        adjacency_list[neg].emplace(pos, true, i);
    }
}
}

#endif // CELL_ADJACENCY_H
