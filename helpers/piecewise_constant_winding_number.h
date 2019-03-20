#ifndef PIECEWISE_CONSTANT_WINDING_NUMBER_H
#define PIECEWISE_CONSTANT_WINDING_NUMBER_H

#include <Eigen/Dense>
#include <vector>

namespace bo {
template <
        typename DerivedF,
        typename DeriveduE,
        typename uE2EType>
inline bool piecewiseConstantWindingNumber(
        const Eigen::MatrixBase<DerivedF>& F,
        const Eigen::MatrixBase<DeriveduE>& uE,
        const std::vector<std::vector<uE2EType> >& uE2E)
{
    const size_t num_faces = F.rows();
    const size_t num_edges = uE.rows();
    const auto edge_index_to_face_index = [&](const size_t &ei) {
        return ei % num_faces;
    };
    const auto is_consistent = [&](const size_t &fid, const size_t &s, const size_t &d) {
        if ((size_t)F(fid, 0) == s && (size_t)F(fid, 1) == d) return true;
        if ((size_t)F(fid, 1) == s && (size_t)F(fid, 2) == d) return true;
        if ((size_t)F(fid, 2) == s && (size_t)F(fid, 0) == d) return true;

        if ((size_t)F(fid, 0) == d && (size_t)F(fid, 1) == s) return false;
        if ((size_t)F(fid, 1) == d && (size_t)F(fid, 2) == s) return false;
        if ((size_t)F(fid, 2) == d && (size_t)F(fid, 0) == s) return false;
        throw std::runtime_error("invalid face in piecewiseConstantWindingNumber");
    };

    for(size_t i = 0; i < num_edges; i++) {
        const size_t &s = uE(i, 0);
        const size_t &d = uE(i, 1);
        int count = 0;
        for(const auto &ei : uE2E[i]) {
            const size_t fid = edge_index_to_face_index(ei);
            if (is_consistent(fid, s, d))
                count++;
            else
                count--;
        }
        if (count != 0)
            return false;
    }
    return true;
}
}

#endif // PIECEWISE_CONSTANT_WINDING_NUMBER_H
