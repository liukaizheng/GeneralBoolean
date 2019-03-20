#ifndef WRITEPATCHES_H
#define WRITEPATCHES_H

#include <Eigen/Core>
#include <vector>
#include "write_obj.h"
#include "assign_scalar.h"

template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedP>
bool writePatches(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const Eigen::MatrixBase<DerivedP>& P)
{
    Eigen::MatrixXd PV;
    PV.resizeLike(V);
    for(Eigen::Index i = 0; i < V.rows(); i++) {
        for(Eigen::Index j = 0; j < V.cols(); j++) {
            bo::assignScalar(V(i, j), PV(i, j));
        }
    }
    const unsigned num_patches = P.maxCoeff() + 1;
    std::vector<std::vector<std::vector<typename DerivedF::Scalar>>> Fs(num_patches);
    for(Eigen::Index i = 0; i <F.rows(); i++) {
        Fs[P(i)].emplace_back(std::vector<typename DerivedF::Scalar>{F(i, 0), F(i, 1), F(i, 2)});
    }

    for(unsigned i = 0; i <num_patches; i++) {
        DerivedF PF(Fs[i].size(), 3);
        for(Eigen::Index j = 0;  j < PF.rows(); j++)
            for(Eigen::Index k = 0; k < PF.cols(); k++)
                PF(j, k) = Fs[i][j][k];
        std::stringstream ss;
        ss << "p_" << i << ".obj";
        std::string name;
        ss >> name;
        bo::writeOBJ(name, V, PF);
    }
    return true;
}

#endif // WRITEPATCHES_H
