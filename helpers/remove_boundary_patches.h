#ifndef RMOVE_BOUNDARY_PATCHES_H_
#define RMOVE_BOUNDARY_PATCHES_H_

#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <algorithm>

template <
    typename DerivedV,
    typename DerivedF,
    typename DerivedP,
    typename DerivedJ,
    typename uE2EType,
    typename DerivedrV,
    typename DerivedrF,
    typename DerivednP,
    typename DerivedrJ>
auto removeBoundaryPatches(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedP> &P,
        const Eigen::MatrixBase<DerivedJ> &J,
        const std::vector<std::vector<uE2EType>> &uE2E,
        Eigen::PlainObjectBase<DerivedrV> &rV,
        Eigen::PlainObjectBase<DerivedrF> &rF,
        Eigen::PlainObjectBase<DerivednP> &nP,
        Eigen::PlainObjectBase<DerivedrJ> &rJ)
{
    const unsigned num_old_faces = static_cast<unsigned>(F.rows());
    const unsigned num_old_patches = static_cast<unsigned>(P.maxCoeff() + 1);
    std::vector<bool> patch_is_boundary(num_old_patches, false);
    for(unsigned i = 0; i < uE2E.size(); i++) {
        if(uE2E[i].size() == 1) {
            const unsigned f = static_cast<unsigned>(uE2E[i][0] % num_old_faces);
            patch_is_boundary[P(f)] = true;
        }
    }

    std::unordered_map<typename DerivedP::Scalar, typename DerivednP::Scalar> op2p;
    typename DerivednP::Scalar patch_count = 0;
    for(typename DerivedP::Scalar i = 0; i < num_old_patches; i++) {
        if(patch_is_boundary[i]) continue;
        op2p.emplace(i, patch_count++);
    }

    std::vector<std::vector<typename DerivedrF::Scalar>> new_faces;
    std::vector<std::vector<typename DerivedrV::Scalar>> new_vertices;
    std::unordered_map<Eigen::Index, Eigen::Index> ov2v;
    std::vector<Eigen::Index> f2of;
    std::vector<typename DerivedrJ::Scalar> f2oof;
    Eigen::Index v_index = 0;
    for(Eigen::Index i = 0; i < F.rows(); i++) {
        if(patch_is_boundary[P(i)]) continue;
        std::vector<typename DerivedF::Scalar> vs;
        vs.reserve(3);
        for(Eigen::Index j = 0; j < 3; j++) {
            const auto &v = F(i, j);
            const auto it = ov2v.find(v);
            if(it == ov2v.end()) {
                vs.emplace_back(v_index++);
                new_vertices.push_back({V(F(i, j), 0), V(F(i, j), 1), V(F(i, j), 2)});
                ov2v.emplace(v, vs[j]);
            }
            else
                vs.emplace_back(it->second);
        }
        new_faces.emplace_back(vs);
        f2of.emplace_back(i);
        f2oof.emplace_back(J(i));
    }

    rV.resize(v_index, 3);
    for(Eigen::Index i = 0; i < v_index; i++)
        for(Eigen::Index j = 0; j < 3; j++)
            rV(i, j) = new_vertices[i][j];

    rF.resize(static_cast<Eigen::Index>(new_faces.size()), 3);
    nP.resize(rF.rows());
    for(Eigen::Index i = 0; i < rF.rows(); i++) {
        nP(i) = op2p[P(f2of[i])];
        for(Eigen::Index j = 0; j < 3; j++)
            rF(i, j) = new_faces[i][j];
    }


    rJ.resize(f2oof.size());
    std::copy(f2oof.begin(), f2oof.end(), rJ.data());
    return patch_count;
}
    

#endif //RMOVE_BOUNDARY_PATCHES_H_
