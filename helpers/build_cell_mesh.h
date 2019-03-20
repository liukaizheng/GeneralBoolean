#ifndef BUILD_CELL_MESH_H_
#define BUILD_CELL_MESH_H_

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "assign_scalar.h"
namespace bo {
template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename PatchType,
    typename CellType,
    typename Index,
    typename DerivedoV,
    typename DerivedoF>
inline void buildCellMesh(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedP> &patch_cells,
        const std::vector<std::vector<PatchType>> &faces_in_patch,
        const std::vector<std::vector<CellType>> &patch_in_cell,
        const Index &i,
        Eigen::PlainObjectBase<DerivedoV> &oV,
        Eigen::PlainObjectBase<DerivedoF> &oF)
{
    typedef typename DerivedF::Scalar FaceIndex;
    std::vector<std::vector<FaceIndex>> faces;
    std::vector<std::vector<typename DerivedV::Scalar>> vertices;
    std::unordered_map<FaceIndex, FaceIndex> ov2v;
    FaceIndex vi = 0;
    for(const auto &p : patch_in_cell[i]) {
        bool reverse = false;
        if(patch_cells(p, 0) == i)
            reverse = true;
        else {
            if(patch_cells(p, 1) != i) throw std::runtime_error("None of the two cells of the patch coicides with the given cell");
        }
        for(const auto &f : faces_in_patch[p]) {
            std::vector<FaceIndex> face;
            if(reverse) face = {F(f, 1), F(f, 0), F(f, 2)};
            else face = {F(f, 0), F(f, 1), F(f, 2)};
            for(auto &v : face) {
                auto it = ov2v.find(v);
                if(it == ov2v.end()) {
                    vertices.push_back({V(v, 0), V(v, 1), V(v, 2)});
                    ov2v.emplace(v, vi);
                    v = vi++;
                }
                else 
                    v = it->second;
            }
            faces.emplace_back(face);
        }
    }

    //std::cout << "face: " << faces.size() << "\n";
    //std::cout << "vertices: " << vertices.size() << "\n";
    oV.resize(static_cast<Eigen::Index>(vi), 3);
    for(Eigen::Index j = 0; j < oV.rows(); j++)
        for(Eigen::Index k = 0; k < 3; k++)
            assignScalar(vertices[j][k], oV(j, k));

    oF.resize(static_cast<Eigen::Index>(faces.size()), 3);
    for(Eigen::Index j = 0; j < oF.rows(); j++)
        for(Eigen::Index k = 0; k < 3; k++)
            oF(j, k) = faces[j][k];
}

}

#endif //BUILD_CELL_MESH_H_
