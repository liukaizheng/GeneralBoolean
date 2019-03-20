#ifndef VERTEX_FACE_ADJACENCY_H
#define VERTEX_FACE_ADJACENCY_H

#include <vector>
#include <Eigen/Dense>

namespace bo {
template <typename SizeT, typename DerivedF, typename Index>
inline void vertexFaceAdjacency(
        const SizeT &n,
        const Eigen::MatrixBase<DerivedF> &F,
        std::vector<std::vector<Index>> &VF,
        std::vector<std::vector<Index>> &VFi)
{
    VF.resize(n);
    VFi.resize(n);

    for(Eigen::Index i = 0; i < F.rows(); i++)
        for(Eigen::Index j = 0; j < F.cols(); j++) {
            VF[F(i,j)].emplace_back(i);
            VFi[F(i,j)].emplace_back(j);
        }
}
}

#endif // VERTEX_FACE_ADJACENCY_H
