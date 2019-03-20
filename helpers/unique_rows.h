#ifndef UNIQUEROWS_H
#define UNIQUEROWS_H

#include "sort_rows.h"
#include <vector>

namespace bo {

/**
 * A: 输入矩阵
 * C: 唯一矩阵
 * IA: 使得C = A(IA)
 * IC: 使得A = C(IC)
 */
template <typename DerivedV, typename DerivedI>
inline void uniqueRows(
        const Eigen::MatrixBase<DerivedV> &A,
        Eigen::PlainObjectBase<DerivedV> &C,
        Eigen::PlainObjectBase<DerivedI> &IA,
        Eigen::PlainObjectBase<DerivedI> &IC)
{
    typedef size_t Index;
    const Index &m = A.rows();
    DerivedV sortA;
    DerivedI IM;
    sortRows(A, sortA, true, IM);
    std::vector<Index> vIC(m);
    for(Index i = 0;i < m; i++)
        vIC[i] = i;
    auto equalIndex = [&](const Index &i, const Index &j) ->bool {
        for(Index c = 0; c < A.cols(); c++) {
            if(sortA(i,c) != sortA(j,c)) return false;
        }
        return true;
    };
    vIC.erase(std::unique(
                  vIC.begin(), vIC.end(), equalIndex), vIC.end());
    IA.resize(m);
    IC.resize(vIC.size());
    Index count = 0;
    IC(0) = IM(0);
    for(Index i = 0; i < m; i++) {
        if(count < vIC.size() - 1 && i == vIC[count+1]) {
            IA(IM(i)) = ++count;
            IC(count) = IM(i);
        }
        else
            IA(IM(i)) = count;
    }
    C.resize(vIC.size(), A.cols());
    for(Index i = 0; i < vIC.size(); i++) {
        C.row(i) = A.row(IC(i));
    }
}
}

#endif // UNIQUEROWS_H
