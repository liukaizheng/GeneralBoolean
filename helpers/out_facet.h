#ifndef OUT_FACET_H
#define OUT_FACET_H

#include "out_edge.h"
#include "order_facets_around_edge.h"

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename IndexType>
inline void outFacet(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedI> &I,
        IndexType &f,
        bool flipped)
{
    typedef typename DerivedF::Scalar Index;
    Index s, d;
    Eigen::Matrix<Index, Eigen::Dynamic, 1> incident_faces;
    outerEdge(V, F, I, s, d, incident_faces);
    assert(incident_faces.size() > 0);

    auto convert_to_signed_index = [&](const Eigen::Index &fid) ->int {
        if ((F(fid, 0) == s && F(fid, 1) == d) ||
                (F(fid, 1) == s && F(fid, 2) == d) ||
                (F(fid, 2) == s && F(fid, 0) == d) ) {
            return (fid+1) * -1;
        } else {
            return fid + 1;
        }
    };

    auto signed_index_to_index = [&](int signed_id) {
        return abs(signed_id) - 1;
    };

    std::vector<int> adj_faces(incident_faces.size());
    std::transform(incident_faces.data(), incident_faces.data() + incident_faces.size(), adj_faces.begin(), convert_to_signed_index);
    Eigen::VectorXi order;
    DerivedV pivot = V.row(s);
    pivot(0) += 1.0;
    orderFacetAroundEdge(V, F, s, d, adj_faces, pivot, order);
    f = signed_index_to_index(adj_faces[order[0]]);
    flipped = adj_faces[order[0]] > 0;
}
}
#endif // OUT_FACET_H
