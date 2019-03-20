#ifndef ORDER_FACETS_AROUND_EDGE_H
#define ORDER_FACETS_AROUND_EDGE_H

#include <Eigen/Dense>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedO>
inline void orderFacetAroundEdge(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const typename DerivedF::Scalar &s,
        const typename DerivedF::Scalar &d,
        const std::vector<int> &adj_faces,
        Eigen::PlainObjectBase<DerivedO> &order)
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
    typedef Kernel::Point_3 Point_3;
    typedef Kernel::Plane_3 Plane_3;
    typedef typename DerivedF::Scalar Index;

    auto get_face_index = [&](const int & signed_fi) ->Index {return abs(signed_fi) - 1;};

    auto get_opposite_vertex = [&](const Index &fid) {
        if(F(fid, 0) != s && F(fid, 0) != d) return F(fid, 0);
        if(F(fid, 1) != s && F(fid, 1) != d) return F(fid, 1);
        if(F(fid, 2) != s && F(fid, 2) != d) return F(fid, 2);
        assert(false);
        return Index(3);
    };

    if(adj_faces.size() == 0) {
        order.resize(0);
        return;
    }
    else if(adj_faces.size() == 1) {
        order.resize(1);
        order[0] = 0;
        return;
    }
    else if(adj_faces.size() == 2) {
        order.resize(2);
        const Index &o1 = get_opposite_vertex(get_face_index(adj_faces[0]));
        const Index &o2 = get_opposite_vertex(get_face_index(adj_faces[1]));
        const Point_3 ps(V(s,0), V(s,1), V(s,2));
        const Point_3 pd(V(d,0), V(d,1), V(d,2));
        const Point_3 p1(V(o1,0), V(o1,1), V(o1,2));
        const Point_3 p2(V(o2,0), V(o2,1), V(o2,2));
        switch (CGAL::orientation(ps, pd, p1, p2)) {
        case CGAL::POSITIVE:
            order(0) = 1;
            order(1) = 0;
            break;
        case CGAL::NEGATIVE:
            order(0) = 0;
            order(1) = 1;
            break;
        case CGAL::COPLANAR:
            switch (CGAL::coplanar_orientation(ps, pd, p1, p2)) {
            case CGAL::POSITIVE:
                if(adj_faces[0] < adj_faces[1]) {
                    order(0) = 0;
                    order(1) = 1;
                }
                else {
                    order(0) = 1;
                    order(1) = 0;
                }
                break;
            case CGAL::NEGATIVE:
                order(0) = 0;
                order(1) = 1;
                break;
            case CGAL::COLLINEAR:
                assert(false && "degenerated triangle detected.\n");
                break;
            default:
                assert(false && "order triangles error\n");
                break;
            }
            break;
        default:
            assert(false && "order triangles error\n");
            break;
        }
        return;
    }

    const Index &num_adj_faces = adj_faces.size();
    const Index o = get_opposite_vertex(get_face_index(adj_faces[0]));
    const Point_3 p_s(V(s,0), V(s,1), V(s,2));
    const Point_3 p_d(V(d,0), V(d,1), V(d,2));
    const Point_3 p_o(V(o,0), V(o,1), V(o,2));
    const Plane_3 seperator(p_s, p_d, p_o);
    if(seperator.is_degenerate())
        assert(false && "degenerated triangle detected");

    std::vector<Point_3> opposite_vertices;
    for(Index i = 0; i < num_adj_faces; i++) {
        const Index o = get_opposite_vertex(get_face_index(adj_faces[i]));
        opposite_vertices.emplace_back(V(o,0), V(o,1), V(o,2));
    }


    std::vector<int> positive_side;
    std::vector<int> negative_side;
    std::vector<int> tie_positive_oriented;
    std::vector<int> tie_negative_oriented;

    std::vector<Index> positive_side_index;
    std::vector<Index> negative_side_index;
    std::vector<Index> tie_positive_oriented_index;
    std::vector<Index> tie_negative_oriented_index;

    for(Index i = 0; i < num_adj_faces; i++) {
        const Point_3 &pa = opposite_vertices[i];
        auto orientation = seperator.oriented_side(pa);
        const Index &f = adj_faces[i];
        switch (orientation) {
        case CGAL::ON_POSITIVE_SIDE:
            positive_side.emplace_back(f);
            positive_side_index.emplace_back(i);
            break;
        case CGAL::ON_NEGATIVE_SIDE:
            negative_side.emplace_back(f);
            negative_side_index.emplace_back(i);
            break;
        case CGAL::ON_ORIENTED_BOUNDARY:
        {
            auto inplane_orientation = CGAL::coplanar_orientation(p_s, p_d, p_o, pa);
            switch (inplane_orientation) {
            case CGAL::POSITIVE:
                tie_positive_oriented.emplace_back(f);
                tie_positive_oriented_index.emplace_back(i);
                break;
            case CGAL::NEGATIVE:
                tie_negative_oriented.emplace_back(f);
                tie_negative_oriented_index.emplace_back(i);
                break;
            case CGAL::COLLINEAR:
                assert(false && "degenerated triangles detected");
                break;
            default:
                assert(false && "order error");
                break;
            }
            break;
        }
        default:
            assert(false && "order error");
            break;
        }
    }

    auto index_sort = [](const std::vector<int> &data) {
        const int len = data.size();
        std::vector<Index> order(len);
        for(int i = 0; i < len; i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](const Index &a, const Index &b) {return data[a] < data[b];});
        return order;
    };
    DerivedO positive_order, negative_order;
    orderFacetAroundEdge(V, F, s, d, positive_side, positive_order);
    orderFacetAroundEdge(V, F, s, d, negative_side, negative_order);
    std::vector<Index> tie_positive_order = index_sort(tie_positive_oriented);
    std::vector<Index> tie_negative_order = index_sort(tie_negative_oriented);

    const Index tie_positive_size = tie_positive_oriented.size();
    const Index tie_negative_size = tie_negative_oriented.size();
    const Index positive_size = positive_order.size();
    const Index negative_size = negative_order.size();

    assert(tie_positive_size + tie_negative_size + positive_size + negative_size == num_adj_faces);

    order.resize(num_adj_faces);

    Index count = 0;
    for(Index i = 0; i < tie_positive_size; i++)
        order(count++) = tie_positive_oriented_index[tie_positive_order[i]];
    for(Index i = 0; i < negative_size; i++)
        order(count++) = negative_side_index[negative_order(i)];
    for(Index i = 0; i < tie_negative_size; i++)
        order(count++) = tie_negative_oriented_index[tie_negative_order[i]];
    for(Index i = 0; i < positive_size; i++)
        order(count++) = positive_side_index[positive_order(i)];

    count = 0;
    for(Index i = 0; i < num_adj_faces - 1; i++) {
        const Point_3 &p_a = opposite_vertices[order(i)];
        const Point_3 &p_b = opposite_vertices[order((i+1)%num_adj_faces)];
        auto orientation = CGAL::orientation(p_s, p_d, p_a, p_b);
        if(orientation == CGAL::POSITIVE) {
            count = (i + 1) % num_adj_faces;
            break;
        }
        else if(orientation == CGAL::COPLANAR &&
                Plane_3(p_s, p_d, p_a).orthogonal_direction() !=
                Plane_3(p_s, p_d, p_b).orthogonal_direction()) {
            count = (i + 1) % num_adj_faces;
            break;
        }
    }
    if(count) {
        DerivedO circular_order = order;
        for(Index i = 0; i < num_adj_faces; i++)
            order(i) = circular_order((i+count) % num_adj_faces);
    }
}

template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename DerivedP,
        typename Index>
inline void orderFacetAroundEdge(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Index &s,
        const Index &d,
        const std::vector<int> &adj_faces,
        const Eigen::MatrixBase<DerivedP> &pivot_point,
        Eigen::PlainObjectBase<DerivedI> &order)
{
    assert(V.cols() == 3);
    assert(F.cols() == 3);
    assert(pivot_point.cols() == 3);
    auto signed_index_to_index = [&](int signed_idx)
    {
        return abs(signed_idx) -1;
    };
    auto get_opposite_vertex_index = [&](const typename DerivedF::Scalar &fid) -> typename DerivedF::Scalar
    {
        typedef typename DerivedF::Scalar IndexType;
        if (F(fid, 0) != (IndexType)s && F(fid, 0) != (IndexType)d) return F(fid, 0);
        if (F(fid, 1) != (IndexType)s && F(fid, 1) != (IndexType)d) return F(fid, 1);
        if (F(fid, 2) != (IndexType)s && F(fid, 2) != (IndexType)d) return F(fid, 2);
        assert(false);
        return 3;
    };

    {
        typedef CGAL::Exact_predicates_exact_constructions_kernel K;
        K::Point_3 ps(V(s,0), V(s,1), V(s,2));
        K::Point_3 pd(V(d,0), V(d,1), V(d,2));
        K::Point_3 pp(pivot_point(0,0), pivot_point(0,1), pivot_point(0,2));
        if (CGAL::collinear(ps, pd, pp)) {
            throw std::runtime_error(
                        "Pivot point is collinear with the outer edge!");
        }
    }

    const size_t &N = adj_faces.size();
    const size_t num_faces = N + 1;

    DerivedV vertices(num_faces + 2, 3);
    for(size_t i = 0; i < N; i++) {
        const auto &fid = signed_index_to_index(adj_faces[i]);
        vertices.row(i) = V.row(get_opposite_vertex_index(fid));
    }
    vertices.row(N) = pivot_point;
    vertices.row(N + 1) = V.row(s);
    vertices.row(N + 2) = V.row(d);

    DerivedF faces(num_faces, 3);
    std::vector<int> adj_faces_with_pivot(num_faces);
    for(size_t i = 0; i < num_faces; i++) {
        if(i == N || adj_faces[i] < 0) {
            faces(i, 0) = N + 1;
            faces(i, 1) = N + 2;
            adj_faces_with_pivot[i] = (i + 1) * -1;
        }
        else {
            faces(i, 0) = N + 2;
            faces(i, 1) = N + 1;
            adj_faces_with_pivot[i] = i + 1;
        }
        faces(i, 2) = i;
    }

    DerivedI order_with_pivot;
    orderFacetAroundEdge(vertices, faces, N + 1, N + 2, adj_faces_with_pivot, order_with_pivot);

    Eigen::Index pivot_index = num_faces;
    for(Eigen::Index i = 0; i < num_faces; i++) {
        if(order_with_pivot[i] == N) {
            pivot_index = i;
            break;
        }
    }

    order.resize(N);
    for(Eigen::Index i = 0; i < N; i++) {
        order(i) = order_with_pivot((i + 1 + pivot_index) % num_faces);
    }
}
}

#endif // ORDER_FACETS_AROUND_EDGE_H
