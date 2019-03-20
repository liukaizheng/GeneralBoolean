#ifndef SORT_ROWS_H
#define SORT_ROWS_H

#include <Eigen/Dense>
#include <algorithm>


namespace bo {
template <typename DerivedV, typename DerivedI>
inline void sortRows(
        const Eigen::MatrixBase<DerivedV> &A,
        Eigen::PlainObjectBase<DerivedV> &C,
        bool is_ascending,
        Eigen::PlainObjectBase<DerivedI> &IC)
{
    typedef Eigen::Index Index;
    const Index &m = A.rows();
    IC.resize(m);
    for(Index i = 0; i < m; i++)
        IC(i) = i;
    if(is_ascending) {
        auto ascending = [&](const Index &i, const Index &j) ->bool {
            for(Index c = 0; c < A.cols(); c++) {
                if(A(i, c) < A(j, c)) return true;
                else if(A(i, c) > A(j, c)) return false;
            }
            return false;
        };
        std::sort(IC.data(), IC.data() + m, ascending);
    }
    else {
        auto descending = [&](const Index &i, const Index &j) ->bool {
            for(Index c = 0; c < A.cols(); c++) {
                if(A(i, c) < A(j, c)) return false;
                else if(A(i, c) > A(j, c)) return true;
            }
            return false;
        };
        std::sort(IC.data(), IC.data() + m, descending);
    }

    C.resize(m, A.cols());
    for(Index i = 0; i < m; i++)
        for(Index j = 0; j < A.cols(); j++)
            C(i, j) = A(IC(i), j);
}
}

#endif // SORT_ROWS_H
