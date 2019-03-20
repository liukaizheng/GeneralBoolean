#ifndef REMOVE_CLOSED_POINTS_H
#define REMOVE_CLOSED_POINTS_H

#include <Eigen/Core>
#include <vector>
#include <unordered_map>
#include <limits>
#include <queue>
#include <type_traits>

template <
        typename DerivedVi,
        typename DerivedFi,
        typename DerivedVo,
        typename DerivedFo,
        typename ScalarType>
void removeClosedPoints(
        const Eigen::MatrixBase<DerivedVi> &Vi,
        const Eigen::MatrixBase<DerivedFi> &Fi,
        ScalarType supreme,
        Eigen::PlainObjectBase<DerivedVo> &Vo,
        Eigen::PlainObjectBase<DerivedFo> &Fo)
{
    typedef typename DerivedFi::Scalar Index;
    std::vector<std::vector<Index>> Fs;
    std::unordered_map<Index, std::vector<Index>> closed_points;
    for(Index i = 0; i < Fi.rows(); i++) {
        bool is_remove = false;
        for(Index j = 0; j < 3; j++) {
            Index vnext = Fi(i, (j+1)%3);
            auto dist = /*static_cast<Eigen::Matrix<typename DerivedVi::Scalar, Eigen::Dynamic, 1>>*/(Vi.row(Fi(i,j)) - Vi.row(vnext)).cwiseAbs().sum();
            if(dist < supreme) {
                is_remove = true;
                auto a_iter = closed_points.find(Fi(i, j));
                if(a_iter == closed_points.end())
                    closed_points.emplace(Fi(i, j), std::vector<Index>{});
                auto b_iter = closed_points.find(vnext);
                if(b_iter == closed_points.end())
                    closed_points.emplace(vnext, std::vector<Index>{});
                closed_points[Fi(i, j)].emplace_back(vnext);
                closed_points[vnext].emplace_back(Fi(i, j));
            }
        }
        if(!is_remove)
            Fs.push_back({Fi(i, 0), Fi(i, 1), Fi(i, 2)});
    }

    Index vcount = 0;
    std::vector<std::remove_reference_t<std::remove_const_t<decltype(Vi.row(0))>>> Vs;
    Index INVALID = std::numeric_limits<Index>::max();
    std::vector<Index>vmap(Vi.rows(), INVALID);
    for(Eigen::Index i = 0; i < Vi.rows(); i++) {
        if(vmap[i] != INVALID) continue;
        vmap[i] = vcount;
        Vs.emplace_back(Vi.row(i));
        std::queue<Index> Q;
        Q.emplace(i);
        while(!Q.empty()) {
            Index cur = Q.front();
            Q.pop();
            auto closed_iter = closed_points.find(cur);
            if(closed_iter != closed_points.end()) {
                for(const auto &c : closed_iter->second) {
                    if(vmap[c] != INVALID) continue;
                    vmap[c] = vcount;
                    Q.emplace(c);
                }
            }
        }
        vcount++;
    }

    Vo.resize(Vs.size(), 3);
    for(Eigen::Index i = 0; i < Vo.rows(); i++)
        Vo.row(i) = Vs[i];

    Fo.resize(Fs.size(), 3);
    for(Eigen::Index i = 0; i < Fo.rows(); i++)
        for(Eigen::Index j = 0; j < 3; j++)
            Fo(i, j) = vmap[Fs[i][j]];

}

#endif // REMOVE_CLOSED_POINTS_H
