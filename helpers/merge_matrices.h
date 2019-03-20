#ifndef MERGE_MATRICES_H
#define MERGE_MATRICES_H

#include <Eigen/Dense>

namespace bo {
template <
        typename DerivedV1,
        typename DerivedF1,
        typename DerivedV2,
        typename DerivedF2,
        typename DerivedVV,
        typename DerivedFF>
inline void mergeMatrices(
        const Eigen::MatrixBase<DerivedV1> &V1,
        const Eigen::MatrixBase<DerivedF1> &F1,
        const Eigen::MatrixBase<DerivedV2> &V2,
        const Eigen::MatrixBase<DerivedF2> &F2,
        Eigen::PlainObjectBase<DerivedVV> &VV,
        Eigen::PlainObjectBase<DerivedFF> &FF)
{
    VV.resize(V1.rows() + V2.rows(), V1.cols());
    VV.block(0, 0, V1.rows(), V1.cols()) = V1;
    VV.block(V1.rows(), 0, V2.rows(), V2.cols()) = V2;

    FF.resize(F1.rows() + F2.rows(), F1.cols());
    FF.block(0, 0, F1.rows(), F1.cols()) = F1;
    FF.block(F1.rows(), 0, F2.rows(), F2.cols()) = (F2.array() + V1.rows()).matrix();
}


}

#endif // MERGE_MATRICES_H
