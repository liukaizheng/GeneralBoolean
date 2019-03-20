#ifndef PROPAGATE_WINDING_NUMBERS_H
#define PROPAGATE_WINDING_NUMBERS_H

//#include "piecewise_constant_winding_number.h"
#include "cell_adjacency.h"
#include <queue>
#include <unordered_map>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DeriveduE,
        typename uE2EType,
        typename DerivedP,
        typename DerivedC,
        typename DerivedL,
        typename DerivedpW,
        typename DerivedW>
inline bool progateWindingNumbers(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DeriveduE> &uE,
        const std::vector<std::vector<uE2EType>> &uE2E,
        const size_t &num_patches,
        const Eigen::MatrixBase<DerivedP> &P,
        const size_t &num_cells,
        const Eigen::MatrixBase<DerivedC> &per_patch_cells,
        const Eigen::MatrixBase<DerivedL> &labels,
        Eigen::PlainObjectBase<DerivedpW> &pW,
        Eigen::PlainObjectBase<DerivedW> &W)
{
    //if(!piecewiseConstantWindingNumber(F, uE, uE2E)) {
        //throw std::runtime_error("input mesh is not PWN");
        //return false;
    //}

    const size_t &num_faces = F.rows();
    typedef std::tuple<typename DerivedC::Scalar, bool, size_t> CellConnection;
    std::vector<std::set<CellConnection>> cell_adj;
    cellAdjacency(per_patch_cells, num_cells, cell_adj);

    DerivedP patch_labels(num_patches);
    const auto &INVALID = std::numeric_limits<int>::max();
    patch_labels.setConstant(INVALID);
    for(size_t i = 0; i < num_faces; i++) {
        if(patch_labels(P(i)) == INVALID) patch_labels(P(i)) = labels(i);
        else assert(patch_labels(P(i)) == labels(i));
    }
    assert((patch_labels.array() != INVALID).all());
    const size_t num_labels = patch_labels.maxCoeff() + 1;

    pW.resize(num_cells, num_labels);
    pW.setConstant(INVALID);
    pW.row(0).setZero();
    std::queue<size_t> Q;
    Q.push(0);
    while(!Q.empty()) {
        const size_t curr_cell = Q.front();
        Q.pop();
        for(const auto &neighbor : cell_adj[curr_cell]) {
            size_t neighbor_cell, patch_idx;
            bool direction;
            std::tie(neighbor_cell, direction, patch_idx) = neighbor;
            if((pW.row(neighbor_cell).array() == INVALID).any()) {
                pW.row(neighbor_cell) = pW.row(curr_cell);
                pW(neighbor_cell, patch_labels(patch_idx)) += direction ? -1 : 1;
                Q.emplace(neighbor_cell);
            }
        }
    }

    W.resize(num_faces, num_labels * 2);
    for(size_t i = 0; i <num_faces; i++) {
        const size_t &patch = P(i);
        const size_t &positive_cell = per_patch_cells(patch, 0);
        const size_t &negative_cell = per_patch_cells(patch, 1);
        for(size_t j = 0; j < num_labels; j++) {
            W(i, 2*j) = pW(positive_cell, j);
            W(i, 2*j+1) = pW(negative_cell, j);
        }
    }

    return true;
}

template <
        typename LabelType,
        typename PatchType,
        typename DerivedF,
        typename DerivedP,
        typename DerivedL,
        typename DerivedC,
        typename CountCellType,
        typename DerivedpW,
        typename DerivedW>

inline bool progateWindingNumbers(
        const std::vector<std::unordered_map<LabelType, std::unordered_map<PatchType, bool>>> &c2l2pd,
        const Eigen::MatrixBase<DerivedF>& F,
        const Eigen::MatrixBase<DerivedP>& P,
        const Eigen::MatrixBase<DerivedL>& patch_labels,
        const Eigen::MatrixBase<DerivedC>& per_patch_cells,
        const CountCellType& num_cells,
        Eigen::PlainObjectBase<DerivedpW>& pW,
        Eigen::PlainObjectBase<DerivedW>& W)
{
    const auto num_faces = F.rows();
    typedef std::tuple<typename DerivedC::Scalar, bool, size_t> CellConnection;
    std::vector<std::set<CellConnection>> cell_adj;
    cellAdjacency(per_patch_cells, num_cells, cell_adj);

    const auto num_labels = patch_labels.maxCoeff() + 1;

    pW.resize(num_cells, num_labels);
    const auto INVALID = std::numeric_limits<int>::max();
    pW.setConstant(INVALID);
    pW.row(0).setConstant(0);

    std::queue<typename DerivedC::Scalar> Q;
    Q.emplace(0);
    while(!Q.empty()) {
        const auto cc = Q.front();
        Q.pop();
        for(const auto& neighbor: cell_adj[cc]) {
            typename DerivedC::Scalar nc;
            size_t p;
            bool direction;
            std::tie(nc, direction, p) = neighbor;
            if((pW.row(nc).array() == INVALID).any()) {
                const auto& cl2pd = c2l2pd[cc];
                const auto& nl2pd = c2l2pd[nc];
                std::unordered_map<LabelType, int> num_non_flippings;
                for(const auto& it: cl2pd) {
                    const auto l = it.first;           //label
                    const auto &pm = it.second;        //patch map
                    const auto pit = pm.find(p);        //patch iterator
                    if(pit != pm.end()) {
                        const auto lit = num_non_flippings.find(l);     //label iterator
                        int inc = pit->second ? 1 : -1;
                        if(lit == num_non_flippings.end())
                            num_non_flippings.emplace(l, std::move(inc));
                        else
                            num_non_flippings[l] += inc;
                    }
                }
                for(const auto& it: nl2pd) {
                    const auto l = it.first;           //label
                    const auto &pm = it.second;        //patch map
                    const auto pit = pm.find(p);        //patch iterator
                    if(pit != pm.end()) {
                        const auto lit = num_non_flippings.find(l);     //label iterator
                        int inc = pit->second ? 1 : -1;
                        if(lit == num_non_flippings.end())
                            num_non_flippings.emplace(l, std::move(inc));
                        else
                            num_non_flippings[l] += inc;
                    }
                }

                for(LabelType i = 0; i <num_labels; i++) {
                    const auto it = num_non_flippings.find(i);
                    if(it == num_non_flippings.end())
                        pW(nc, i) = pW(cc, i);
                    else {
                        pW(nc, i) = pW(cc, i) + (-it->second *(direction ? -1 : 1));
                    }
                }
                Q.emplace(nc);
            }
        }
    }

    W.resize(num_faces, num_labels * 2);
    for(size_t i = 0; i <num_faces; i++) {
        const size_t &patch = P(i);
        const size_t &positive_cell = per_patch_cells(patch, 0);
        const size_t &negative_cell = per_patch_cells(patch, 1);
        for(size_t j = 0; j < num_labels; j++) {
            W(i, 2*j) = pW(positive_cell, j);
            W(i, 2*j+1) = pW(negative_cell, j);
        }
    }

    return true;
}
}

#endif // PROPAGATE_WINDING_NUMBERS_H
