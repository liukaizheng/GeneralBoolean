#ifndef EXTRACT_PATCH_ADJACENCY_H
#define EXTRACT_PATCH_ADJACENCY_H

#include <Eigen/Core>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

template <
        typename DerivedP,
        typename DerivedE,
        typename uE2EType,
        typename CountPatchType,
        typename PatchType>

inline void extractPatchAdjaceny(
        const Eigen::MatrixBase<DerivedP>& P,
        const Eigen::MatrixBase<DerivedE>& E,
        const std::vector<std::vector<uE2EType>>& uE2E,
        const CountPatchType& num_patches,
        std::vector<std::unordered_set<PatchType>>& P2P,
        std::vector<std::unordered_map<PatchType, bool>> &P2P_D)
{
    const auto num_faces = E.rows() / 3;
    std::cout << "num_faces: " << num_faces << std::endl;
    P2P.resize(num_patches);
    P2P_D.resize(num_patches);
    for(const auto&es: uE2E) {
        if(es.size() > 1) {
            std::vector<typename DerivedP::Scalar> ps;
            std::vector<typename DerivedE::Scalar> first_verts;
            for(unsigned i = 0; i < es.size(); i++) {
                const auto fi = es[i] % num_faces;
                const auto p = P(fi);
                const auto it = std::find(ps.begin(), ps.end(), p);
                if(it == ps.end()) {
                    ps.emplace_back(p);
                    first_verts.emplace_back(E(es[i], 0));
                }
            }
            if(ps.size() > 1) {
                for(unsigned i = 0; i < ps.size(); i++) {
                    const auto p = ps[i];
                    for(unsigned j = 1; j < ps.size(); j++) {
                        const auto np = ps[j];
                        bool consisent = first_verts[i] != first_verts[j] ? true : false;
                        P2P[p].emplace(np);
                        P2P[np].emplace(p);
                        P2P_D[p].emplace(np, consisent);
                        P2P_D[np].emplace(p, consisent);
                    }
                }
            }
        }
    }
//    for(const auto& es: uE2E) {
//        if(es.size() > 2) {
//            for(unsigned i = 0; i < es.size(); i++) {
//                const auto fi = es[i] % num_faces;
//                const auto &ei0 = E(es[i], 0);
//                for(unsigned j = i + 1; j <es.size(); j++) {
//                    const auto fj = es[j] % num_faces;
//                    const auto &ej0 = E(es[j], 0);
//                    const auto &ej1 = E(es[j], 1);
//                    P2P[P(fi)].emplace(P(fj));
//                    P2P[P(fj)].emplace(P(fi));
//                    bool consistent = false;
//                    if(ei0 != ej0) {
//                        assert(ei0 == ej1);
//                        consistent = true;
//                    }
//                    P2P_D[P(fi)].emplace(P(fj), consistent);
//                    P2P_D[P(fj)].emplace(P(fi), consistent);
//                }
//            }
//        }
//    }
}

#endif // EXTRACT_PATCH_ADJACENCY_H
