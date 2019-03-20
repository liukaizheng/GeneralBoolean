#ifndef EXTRACT_COMPONENTS_H
#define EXTRACT_COMPONENTS_H

#include <Eigen/Dense>
#include <vector>
#include <queue>

namespace bo {
template <typename TTIndex, typename DerivedC>
inline size_t extractComponents(
        const std::vector<std::vector<std::vector<TTIndex>>> &TT,
        Eigen::PlainObjectBase<DerivedC> &C)
{
    typedef typename DerivedC::Scalar Index;
    const size_t &m = TT.size();
    C.resize(m);
    C.setConstant(m);
    size_t count = 0;
    for(size_t i = 0; i < m; i++) {
        if(C(i) != m) continue;
        std::queue<size_t> Q;
        Q.push(i);
        C(i) = count;
        while(!Q.empty()) {
            const auto idx = Q.front();
            Q.pop();
            for(size_t j = 0; j < 3; j++) {
                for(size_t k = 0; k < TT[idx][j].size(); k++) {
                    const auto &cur = TT[idx][j][k];
                    if(C(cur) != m) continue;
                    C(cur) = count;
                    Q.push(cur);
                }
            }
        }
        count++;
    }
    return count;
}
}

#endif // EXTRACT_COMPONENTS_H
