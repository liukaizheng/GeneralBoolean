#ifndef REMOVE_UNREFERENCED_VERTICES_H
#define REMOVE_UNREFERENCED_VERTICES_H

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedVV,
        typename DerivedFF>
inline void removeUnreferencedVertices(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        Eigen::PlainObjectBase<DerivedVV> &VV,
        Eigen::PlainObjectBase<DerivedFF> &FF)
{
    const size_t &num_verts = V.rows();
    std::vector<size_t> counts(num_verts, 0);
    const size_t &num_faces = F.rows();
    for(size_t i = 0; i < num_faces; i++)
        for(size_t j = 0; j < F.cols(); j++) {
            counts[F(i,j)]++;
        }
    std::unordered_map<size_t, size_t> old_to_new;
    std::vector<size_t> new_to_old;
    size_t num_new_verts = 0;
    for(size_t i = 0; i < num_verts; i++) {
        if(counts[i]) {
            old_to_new.emplace(i, num_new_verts);
            new_to_old.emplace_back(i);
            num_new_verts++;
        }
    }

    VV.resize(num_new_verts, V.cols());
    for(size_t i = 0; i < num_new_verts; i++)
        VV.row(i) = V.row(new_to_old[i]);

    FF.resizeLike(F);
    for(size_t i = 0; i < num_faces; i++)
        for(size_t j = 0; j < F.cols(); j++)
            FF(i, j) = old_to_new[F(i, j)];
}
}

#endif // REMOVE_UNREFERENCED_VERTICES_H
