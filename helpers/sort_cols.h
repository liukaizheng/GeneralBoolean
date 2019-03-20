#ifndef SORT_COLS_H
#define SORT_COLS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace bo {
template <typename DerivedA>
inline void sortCols(
        const Eigen::MatrixBase<DerivedA> &A,
        const bool &is_ascending,
        Eigen::PlainObjectBase<DerivedA> &C)
{
    typedef typename DerivedA::Scalar Scalar;
    C.resize(A.rows(), A.cols());
    std::vector<Scalar> ss(A.cols());
    for(size_t i = 0; i < A.rows(); i++) {
        for(size_t j = 0; j < A.cols(); j++)
            ss[j] = A(i, j);
        std::sort(ss.begin(), ss.end(), [&](const Scalar a, const Scalar b) {
            if(is_ascending) return a < b;
            else return a > b;
        });
        for(size_t j = 0; j < A.cols(); j++)
            C(i, j) = ss[j];
    }
}
}

#endif // SORT_COLS_H
