#ifndef LISTTOMATRIX_H
#define LISTTOMATRIX_H
#include <vector>
#include <Eigen/Dense>

namespace bo {
template <typename scalar, typename Derived>
inline bool listToMatrix(const std::vector<std::vector<scalar>> &V, Eigen::PlainObjectBase<Derived> &M)
{
    assert(V.size() > 0);
    assert(V[0].size() > 0);
    M.resize(V.size(), V[0].size());
    for(size_t i = 0; i < V.size(); i++)
        for(size_t j = 0; j < V[0].size(); j++)
            M(i, j) = V[i][j];
    return true;
}
}
#endif // LISTTOMATRIX_H
