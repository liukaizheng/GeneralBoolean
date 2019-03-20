#ifndef EXTRACT_MANIFOLD_PATCHES_H
#define EXTRACT_MANIFOLD_PATCHES_H

#include <Eigen/Dense>
#include <vector>
#include <queue>

namespace bo {
template <typename DerivedF, typename DerivedI, typename uE2EType>
inline typename DerivedI::Index extractManifoldPatches(
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedI> &EMAP,
        const std::vector<std::vector<uE2EType>> &uE2E,
        Eigen::PlainObjectBase<DerivedI> &P)
{
    typedef typename DerivedI::Scalar Index;
    const Index &num_faces = F.rows();
    assert(F.cols() == 3);
    auto edge_to_face = [&](const Index &ei) ->Index { return ei % num_faces;};
    auto face_and_corner_to_edge = [&](const Index &fi, const Index &ci) ->Index {return ci * num_faces + fi;};
    auto is_manifold_edge = [&](const Index &fi, const Index &ci) {
        const auto &ei = face_and_corner_to_edge(fi, ci);
        return uE2E[EMAP(ei)].size() == 2;
    };
    auto get_adj_face = [&](const Index &fi, const Index &ci) ->Index {
        const auto &ei = face_and_corner_to_edge(fi, ci);
        const auto &adj_faces = uE2E[EMAP(ei)];
        assert(adj_faces.size() == 2);
        if(adj_faces[0] == ei) return edge_to_face(adj_faces[1]);
        else return edge_to_face(adj_faces[0]);
    };

    const Index INVALID = num_faces;
    P.resize(num_faces);
    P.setConstant(INVALID);
    Index num_patches = 0;
    for(Index i = 0; i < num_faces; i++) {
        if(P(i) != INVALID) continue;

        std::queue<Index> Q;
        Q.push(i);
        P(i) = num_patches;
        while(!Q.empty()) {
            const Index cur = Q.front();
            Q.pop();
            for(Index i = 0; i < 3; i++) {
                if(is_manifold_edge(cur, i)) {
                    const Index adj_face = get_adj_face(cur, i);
                    if(P(adj_face) != INVALID) continue;
                    Q.push(adj_face);
                    P(adj_face) = num_patches;
                }
            }
        }
        num_patches++;
    }
    assert((P.array() != INVALID).all());
    return num_patches;
}
}

#endif // EXTRACT_MANIFOLD_PATCHES_H
