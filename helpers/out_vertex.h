#ifndef OUT_VERTEX_H
#define OUT_VERTEX_H

#include <vector>
#include <Eigen/Dense>
#include <algorithm>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename IndexType,
        typename DerivedA>
inline void outerVertex(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedI> &I,
        IndexType &v_index,
        Eigen::PlainObjectBase<DerivedA> &A)
{
    typedef typename DerivedF::Scalar Index;
    const Index INVALID = std::numeric_limits<Index>::max();
    const Eigen::Index num_selected_faces = I.size();
    std::vector<typename DerivedA::Scalar> candidate_faces;
    Index out_vid = INVALID;
    typename DerivedV::Scalar out_val;
    for(Eigen::Index i = 0; i < num_selected_faces; i++) {
        const typename DerivedI::Scalar &f = I(i);
        for(Eigen::Index i = 0; i < 3; i++) {
            const Index &v = F(f, i);
            const auto &vx = V(v, 0);
            if(out_vid == INVALID || vx > out_val) {
                out_val = vx;
                out_vid = v;
                candidate_faces = {f};
            }
            else if(v == out_vid) {
                candidate_faces.emplace_back(f);
            }
            else if(vx == out_val) {
                const auto &vy = V(v, 1);
                const auto &vz = V(v, 2);
                const auto &out_y = V(out_vid, 1);
                const auto &out_z = V(out_vid, 2);
                bool replace = (vy > out_y ||
                                ((vy == out_y) && (vz > out_z)));
                if(replace) {
                    out_val = vx;
                    out_vid = v;
                    candidate_faces = {f};
                }
            }
        }
    }

    v_index = out_vid;
    A.resize(candidate_faces.size());
    std::copy(candidate_faces.begin(), candidate_faces.end(), A.data());
}
}

#endif // OUT_VERTEX_H
