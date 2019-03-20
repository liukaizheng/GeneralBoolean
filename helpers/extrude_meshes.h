#include <Eigen/Dense>
#include <vector>
#include <iostream>

template <
    typename DerivedV,
    typename Scalar>
inline void extrudeMeshes(
        const std::vector<DerivedV> &Vs,
        const Scalar &scale,
        std::vector<DerivedV> &oVs)
{
    const auto &num_meshes = Vs.size();
    std::vector<Eigen::Matrix<typename DerivedV::Scalar, 1, 3>> mesh_centers(num_meshes);
    Eigen::Index num_verts = 0;
    for(unsigned i = 0; i < num_meshes; i++) {
        mesh_centers[i] = Vs[i].colwise().sum() / Vs[i].rows();
        num_verts += Vs[i].rows();
    }

    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> center;
    center.setZero();
    for(unsigned i = 0; i < num_meshes; i++) {
        center +=  mesh_centers[i] * (static_cast<Scalar>(Vs[i].rows()) / num_verts);
    }

    oVs.resize(num_meshes);
    for(unsigned i = 0 ; i < num_meshes; i++) {
        auto distance = (mesh_centers[i] - center) * scale;
        oVs[i].resizeLike(Vs[i]);
        for (unsigned j = 0; j < Vs[i].rows(); j++) {
            oVs[i].row(j) = Vs[i].row(j) + distance;
        }
    }
    //std::cout << "mesh_centers[0]  ";
    //std::cout << mesh_centers[0] << "\n";

}
