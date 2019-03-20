#ifndef RESOLVE_DUPLICATE_FACES_H
#define RESOLVE_DUPLICATE_FACES_H

#include "unique_rows.h"

namespace bo {

/**
 * F2 = F1(J)
 */
template <typename DerivedF1, typename DerivedF2, typename DerivedJ>
inline void resolveDuplicateFaces(
        const Eigen::MatrixBase<DerivedF1> &F1,
        Eigen::PlainObjectBase<DerivedF2> &F2,
        Eigen::PlainObjectBase<DerivedJ> &J)
{
    DerivedF1 sortF = F1;
    for(Eigen::Index i = 0; i < F1.rows(); i++) {
        std::vector<typename DerivedF1::Scalar *> f(3);
        for(Eigen::Index j = 0; j < 3; j++)
            f[j] = &sortF(i, j);
        if(*f[0] > *f[1]) std::swap(*f[0], *f[1]);
        if(*f[1] > *f[2]) std::swap(*f[1], *f[2]);
        if(*f[0] > *f[1]) std::swap(*f[0], *f[1]);
    }
    DerivedF1 uF;
    DerivedJ IA, IC;
    uniqueRows(sortF, uF, IC, IA); //uF = sortF(IC)
//    if(uF.rows() == F1.rows()) {
//        F2 = F1;
//        J = DerivedJ::LinSpaced(F1.rows(), 0, F1.rows() - 1);
//        return;
//    }

    const size_t &num_faces = F1.rows();
    const size_t &num_unique_faces = uF.rows();
    Eigen::VectorXi &&counts = Eigen::VectorXi::Zero(num_unique_faces);
    Eigen::VectorXi &&ucounts = Eigen::VectorXi::Zero(num_unique_faces);
    std::vector<std::vector<int>> uF2F(num_unique_faces);
    for(size_t i = 0; i < num_faces; i++) {
        const auto &ui = IC(i);
        const bool consistent =
                (F1(i,0) == uF(ui, 0) && F1(i,1) == uF(ui, 1) && F1(i,2) == uF(ui, 2)) ||
                (F1(i,0) == uF(ui, 1) && F1(i,1) == uF(ui, 2) && F1(i,2) == uF(ui, 0)) ||
                (F1(i,0) == uF(ui, 2) && F1(i,1) == uF(ui, 0) && F1(i,2) == uF(ui, 1));
        uF2F[ui].emplace_back(int(i + 1) * (consistent? 1 : -1));
        ucounts(ui)++;
        counts(ui) += consistent ? 1 : -1;
    }

    std::vector<size_t> kept_faces;
    for (size_t i=0; i<num_unique_faces; i++) {
        if (ucounts[i] == 1) {
            kept_faces.push_back(abs(uF2F[i][0])-1);
            continue;
        }
        if (counts[i] == 1) {
            bool found = false;
            for (const auto &fid : uF2F[i]) {
                if (fid > 0) {
                    kept_faces.push_back(abs(fid)-1);
                    found = true;
                    break;
                }
            }
            assert(found);
        } else if (counts[i] == -1) {
            bool found = false;
            for (const auto &fid : uF2F[i]) {
                if (fid < 0) {
                    kept_faces.push_back(abs(fid)-1);
                    found = true;
                    break;
                }
            }
            assert(found);
        } else {
            assert(counts[i] == 0);
        }
    }
    const size_t &num_kept = kept_faces.size();
    J.resize(num_kept);
    std::copy(kept_faces.begin(), kept_faces.end(), J.data());

    F2.resize(num_kept, F1.cols());
    for(size_t i = 0; i < num_kept; i++)
        F2.row(i) = F1.row(J(i));

}

}

#endif // RESOLVE_DUPLICATE_FACES_H
