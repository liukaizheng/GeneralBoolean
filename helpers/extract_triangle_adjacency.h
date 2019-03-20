#ifndef EXTRACT_TRIANGLE_ADJACENCY_H
#define EXTRACT_TRIANGLE_ADJACENCY_H

#include <Eigen/Dense>
#include <vector>
#include <algorithm>

namespace bo {
template <
        typename DerivedEMAP,
        typename uE2EType>
inline void extractTriangleAdjcency(
        const Eigen::MatrixBase<DerivedEMAP> &EMAP,
        const std::vector<std::vector<uE2EType>> &uE2E,
        std::vector<std::vector<std::vector<size_t>>> &TT)
{
    typedef typename DerivedEMAP::Scalar Index;
    assert(EMAP.size() % 3 == 0);
    const auto num_faces = EMAP.size() / 3;
    TT.resize(num_faces);
    for(Eigen::Index i = 0; i < num_faces; i++) {
        TT[i].resize(3);
        for(Index j = 0; j < 3; j++)
            TT[i][j].reserve(1);
    }
    auto edge_to_face = [&](const size_t &ei) {return ei % num_faces;};
    auto get_edge_index = [&](const size_t &fid, const size_t &idx) {return idx * num_faces + fid;};

    for(Eigen::Index i = 0; i < num_faces; i++) {
        for(size_t j = 0; j < 3; j++) {
            const size_t ei = get_edge_index(i, j);
            const auto &adj_edges = uE2E[EMAP(ei)];
            for(size_t k = 0; k < adj_edges.size(); k++) {
                const size_t fi = edge_to_face(adj_edges[k]);
                if(fi != i) {
                    if(TT[i][j].empty()) TT[i][j].emplace_back(fi);
                    else {
                        const auto &it = std::find(TT[i][j].begin(), TT[i][j].end(), fi);
                        if(it == TT[i][j].end()) TT[i][j].emplace_back(fi);
                    }
                }

            }
        }
    }
}
}

#endif // EXTRACT_TRIANGLE_ADJACENCY_H
