#ifndef RESIZE_MESHES_H
#define RESIZE_MESHES_H

#include <Eigen/Dense>
#include <vector>
#include <type_traits>

template <
        typename DerivedV,
        typename Scalar>
inline void resizeMeshes(
        const Eigen::PlainObjectBase<DerivedV> &V,
        const Scalar &s,
        Eigen::PlainObjectBase<DerivedV> &oV)
{
    oV.resizeLike(V);
    const auto &num_verts = V.rows();
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> center;
    Eigen::Array<typename DerivedV::Scalar, 1, 3> min_row, max_row;
    center.setZero();
    max_row = min_row = V.row(0).array();
    for(Eigen::Index i = 0; i < num_verts; i++) {
        center += V.row(i);
        auto a = V.row(i).array();
        max_row = max_row.max(a);
        min_row = min_row.min(a);
    }
    center /= num_verts;
    auto scale = (max_row - min_row).maxCoeff();
    scale = s / scale;
    for(Eigen::Index i = 0; i < V.rows(); i++) {
        oV.row(i) = (V.row(i).matrix() - center) * scale;
    }
//    std::cout << "max row: " << max_row << "\n";
}

#endif // RESIZE_MESHES_H
