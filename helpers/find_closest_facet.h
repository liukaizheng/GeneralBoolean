#ifndef FIND_CLOSEST_FACET_H
#define FIND_CLOSEST_FACET_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <set>
#include <unordered_map>
#include <algorithm>
#include "order_facets_around_edge.h"

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename DerivedP,
        typename uE2EType,
        typename DerivedEMAP,
        typename Kernel,
        typename VFType,
        typename VFiType>
inline auto findClosetFacet(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedI> &I,
        const Eigen::MatrixBase<DerivedP> &P,
        const std::vector<std::vector<uE2EType>> &uE2E,
        const Eigen::MatrixBase<DerivedEMAP> &EMAP,
        const std::vector<std::vector<VFType>> &VF,
        const std::vector<std::vector<VFiType>> &VFi,
        const CGAL::AABB_tree<
            CGAL::AABB_traits<
                Kernel,
                CGAL::AABB_triangle_primitive<
                    Kernel,
                    typename std::vector<typename Kernel::Triangle_3>::iterator>>> &tree,
        const std::vector<typename Kernel::Triangle_3> &triangles,
        const std::vector<bool> in_I,
        size_t &closest_face,
        bool &is_flip)
{
    if (F.rows() <= 0 || I.rows() <= 0) {
        throw std::runtime_error(
                    "Closest facet cannot be computed on empty mesh.");
    }

    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Segment_3 Segment_3;
    typedef typename Kernel::Triangle_3 Triangle_3;
    typedef typename Kernel::Plane_3 Plane_3;
    typedef typename Kernel::Direction_3 Direction_3;
    typedef typename Kernel::Vector_3 Vector_3;
    typedef typename Kernel::Line_3 Line_3;
    typedef typename Kernel::FT FT;

    enum ElementType { VERTEX, EDGE, FACE };
    auto determine_element_type = [&](const Point_3& p, const size_t fid,
            size_t& element_index) -> ElementType {
        const auto& tri = triangles[fid];
        const Point_3 &p0 = tri[0];
        const Point_3 &p1 = tri[1];
        const Point_3 &p2 = tri[2];

        if (p == p0) { element_index = 0; return VERTEX; }
        if (p == p1) { element_index = 1; return VERTEX; }
        if (p == p2) { element_index = 2; return VERTEX; }
        if (CGAL::collinear(p0, p1, p)) { element_index = 2; return EDGE; }
        if (CGAL::collinear(p1, p2, p)) { element_index = 0; return EDGE; }
        if (CGAL::collinear(p2, p0, p)) { element_index = 1; return EDGE; }

        element_index = 0;
        return FACE;
    };

    auto index_to_signed_index = [&](const uE2EType &index, bool ori) -> int{
        return (index+1) * (ori? 1:-1);
    };

    auto get_orientation = [&](const uE2EType &fid, size_t s, size_t d) -> bool
    {
        const auto& f = F.row(fid);
        if      ((size_t)f[0] == s && (size_t)f[1] == d) return false;
        else if ((size_t)f[1] == s && (size_t)f[2] == d) return false;
        else if ((size_t)f[2] == s && (size_t)f[0] == d) return false;
        else if ((size_t)f[0] == d && (size_t)f[1] == s) return true;
        else if ((size_t)f[1] == d && (size_t)f[2] == s) return true;
        else if ((size_t)f[2] == d && (size_t)f[0] == s) return true;
        else {
            throw std::runtime_error(
                        "Cannot compute orientation due to incorrect connectivity");
            return false;
        }
    };
    auto on_the_positive_side = [&](const uE2EType &fid, const Point_3& p) -> bool
    {
        const auto& f = F.row(fid).eval();
        Point_3 v0(V(f[0], 0), V(f[0], 1), V(f[0], 2));
        Point_3 v1(V(f[1], 0), V(f[1], 1), V(f[1], 2));
        Point_3 v2(V(f[2], 0), V(f[2], 1), V(f[2], 2));
        auto ori = CGAL::orientation(v0, v1, v2, p);
        switch (ori) {
        case CGAL::POSITIVE:
            return true;
        case CGAL::NEGATIVE:
            return false;
        case CGAL::COPLANAR:
            return false;
        default:
            throw std::runtime_error("Unknown CGAL state.");
        }
        return false;
    };

    auto process_edge_case =[&] (
            const Point_3 &query_point,
            const typename DerivedF::Scalar &s,
            const typename DerivedF::Scalar &d,
            size_t preferred_facet,
            bool &orientation) {
        size_t corner_idx = std::numeric_limits<size_t>::max();
        if ((s == F(preferred_facet, 0) && d == F(preferred_facet, 1)) ||
                (s == F(preferred_facet, 1) && d == F(preferred_facet, 0)))
        {
            corner_idx = 2;
        } else if ((s == F(preferred_facet, 0) && d == F(preferred_facet, 2)) ||
                   (s == F(preferred_facet, 2) && d == F(preferred_facet, 0)))
        {
            corner_idx = 1;
        } else if ((s == F(preferred_facet, 1) && d == F(preferred_facet, 2)) ||
                   (s == F(preferred_facet, 2) && d == F(preferred_facet, 1)))
        {
            corner_idx = 0;
        } else
        {
            throw std::runtime_error(
                        "Invalid connectivity, edge does not belong to facet");
        }

        const auto &ueid = EMAP(preferred_facet + corner_idx * F.rows());
        const auto &eids = uE2E[ueid];
        std::vector<uE2EType> intersected_face_indices;

        for(const auto &eid : eids) {
            const auto fid = eid % F.rows();
            if(in_I[fid]) {
                intersected_face_indices.emplace_back(fid);
            }
        }

        const size_t &num_intersected_faces = intersected_face_indices.size();
        std::vector<int> intersected_face_signed_indices(num_intersected_faces);

        std::transform(
                    intersected_face_indices.begin(),
                    intersected_face_indices.end(),
                    intersected_face_signed_indices.begin(),
                    [&](const uE2EType &index) {return index_to_signed_index(index, get_orientation(index, s, d));});

        assert(num_intersected_faces > 0);

        if(num_intersected_faces < 2) {
            const auto &fid = intersected_face_indices[0];
            orientation = on_the_positive_side(fid, query_point);
            return fid;
        }

        DerivedI order;
        orderFacetAroundEdge(V, F, s, d, intersected_face_signed_indices, P, order);

        const size_t first = order[0];
        const size_t last = order[num_intersected_faces - 1];
        if(preferred_facet == intersected_face_indices[first]) {
            orientation = intersected_face_signed_indices[first] > 0;
            return intersected_face_indices[first];
        }
        else if(preferred_facet == intersected_face_indices[last]) {
            orientation = intersected_face_signed_indices[last] < 0;
            return intersected_face_indices[last];
        }
        else {
            orientation = intersected_face_signed_indices[first] > 0;
            return intersected_face_indices[first];
        }
    };

    auto process_vertex_case = [&] (
            const Point_3 &query_point,
            const typename DerivedF::Scalar &s,
            bool &orientation)
    {
        const Point_3 closest_point(V(s,0), V(s,1), V(s,2));
        const auto &all_adj_faces = VF[s];

        std::vector<VFType> adj_faces;
        for(unsigned i = 0; i < all_adj_faces.size(); i++) {
            const auto &fid = all_adj_faces[i];
            if(in_I[fid]) adj_faces.emplace_back(fid);
        }

//        const Direction_3 direction(closest_point, query_point);
        const Plane_3 seperator(query_point, Direction_3(Line_3(closest_point, query_point)));
        const Vector_3 normal = seperator.orthogonal_direction().vector();
        const auto &num_of_adj_faces = adj_faces.size();

        std::vector<FT> angles(num_of_adj_faces);
        std::vector<bool> is_negetive(num_of_adj_faces, false);
        for(unsigned i = 0; i < num_of_adj_faces; i++) {
            const auto &f = all_adj_faces[i];
            const Plane_3 face_plane(
                        Point_3(V(F(f,0), 0), V(F(f,0), 1), V(F(f,0), 2)),
                        Point_3(V(F(f,1), 0), V(F(f,1), 1), V(F(f,1), 2)),
                        Point_3(V(F(f,2), 0), V(F(f,2), 1), V(F(f,2), 2)));
            const Vector_3 face_normal = face_plane.orthogonal_direction().vector();
            angles[i] = normal * face_normal;
            if(CGAL::is_negative(angles[i])) {
                angles[i] = CGAL::abs(angles[i]);
                is_negetive[i] = true;
            }
        }
        std::ptrdiff_t ind = std::max_element(angles.begin(), angles.end()) - angles.begin();
        orientation = is_negetive[ind];
        return adj_faces[ind];
    };

    /*auto process_vertex_case = [&] (
            const Point_3 &query_point,
            const typename DerivedF::Scalar &s,
            const size_t preferred_facet,
            bool &orientation) {
        const Point_3 closest_point(V(s,0), V(s,1), V(s,2));
        std::vector<size_t> adj_faces;
        std::vector<size_t> adj_face_corners;

        const auto &all_adj_faces = VF[s];
        const auto &all_adj_face_corners = VFi[s];
        const size_t &num_all_adj_faces = all_adj_faces.size();
        for(size_t i = 0; i < num_all_adj_faces; i++) {
            const size_t fid = all_adj_faces[i];
            if(in_I[fid]) {
                adj_faces.emplace_back(fid);
                adj_face_corners.emplace_back(all_adj_face_corners[i]);
            }
        }

        const size_t &num_adj_faces = adj_faces.size();
        assert(num_adj_faces > 0);

        std::set<typename DerivedF::Scalar> adj_vertices_set;
        std::unordered_multimap<typename DerivedF::Scalar, size_t> v2f;
        for(size_t i = 0; i < num_adj_faces; i++) {
            const size_t &fid = adj_faces[i];
            const size_t &cid = adj_face_corners[i];
            const auto &f = F.row(fid);
            const auto &next = f[(cid+1) % 3];
            const auto &prev = f[(cid+2) % 3];
            adj_vertices_set.insert({next, prev});
            v2f.insert({{next, fid}, {prev, fid}});
        }

        const size_t &num_adj_vertices = adj_vertices_set.size();
        std::vector<typename DerivedF::Scalar> adj_vertices(num_adj_vertices);
        std::copy(adj_vertices_set.begin(), adj_vertices_set.end(), adj_vertices.begin());

        std::vector<Point_3> adj_points;
        for(const auto &vid : adj_vertices) {
            adj_points.emplace_back(V(vid,0), V(vid,1), V(vid,2));
        }

        auto is_on_exterior = [&](const Plane_3& separator) -> bool{
            size_t positive=0;
            size_t negative=0;
            size_t coplanar=0;
            for (const auto& point : adj_points) {
                switch(separator.oriented_side(point)) {
                case CGAL::ON_POSITIVE_SIDE:
                    positive++;
                    break;
                case CGAL::ON_NEGATIVE_SIDE:
                    negative++;
                    break;
                case CGAL::ON_ORIENTED_BOUNDARY:
                    coplanar++;
                    break;
                default:
                    throw "Unknown plane-point orientation";
                }
            }
            auto query_orientation = separator.oriented_side(query_point);
            if (query_orientation == CGAL::ON_ORIENTED_BOUNDARY &&
                    (positive == 0 && negative == 0)) {
                return true;
            } else {
                bool r = (positive == 0 && query_orientation == CGAL::POSITIVE)
                        || (negative == 0 && query_orientation == CGAL::NEGATIVE);
                return r;
            }
        };

        typename DerivedF::Scalar d = std::numeric_limits<int>::max();
        for(size_t i = 0; i < num_adj_vertices; i++) {
            const auto &vi = adj_vertices[i];
            for(size_t j = i+1; j < num_adj_vertices; j++) {
                Plane_3 seperator(closest_point, adj_points[i], adj_points[j]);
                if(seperator.is_degenerate()) continue;

                if(is_on_exterior(seperator)) {
                    if(!CGAL::collinear(query_point, adj_points[i], closest_point)) {
                        d = vi;
                        break;
                    }
                    else {
                        d = adj_vertices[j];
                        assert(!CGAL::collinear(query_point, adj_points[j], closest_point));
                        break;
                    }
                }
            }
        }

        assert(d != std::numeric_limits<int>::max());

        const auto &itr = v2f.equal_range(d);
        assert(itr.first != itr.second);
        return process_edge_case(query_point, s, d, itr.first->second, orientation);
    };*/

    auto process_face_case = [&] (
            const Point_3 &query_point,
            const size_t &fid,
            bool &orientation) {
        const auto &f = F.row(I(fid));
        return process_edge_case(query_point, f(0), f(1), I(fid), orientation);
    };

    const Point_3 query(P(0), P(1), P(2));
    auto projection = tree.closest_point_and_primitive(query);
    const Point_3 &closest_point = projection.first;
    size_t fid = projection.second - triangles.begin();

    size_t element_index;
    auto element_type = determine_element_type(closest_point, fid, element_index);

    switch (element_type) {
    case VERTEX: {
        const auto &f = F.row(I(fid));
        const auto &s = f[element_index];
        closest_face = process_vertex_case(query, s, is_flip);
    }
        break;
    case EDGE: {
        const auto &f = F.row(I(fid));
        const auto &s = f((element_index+1)%3);
        const auto &d = f((element_index+2)%3);
        closest_face = process_edge_case(query, s, d, I(fid), is_flip);
    }
        break;
    case FACE:
        closest_face = process_face_case(query, fid, is_flip);
        break;
    default:
        assert("unkown element type!");
        break;
    }

    return tree.squared_distance(query, closest_point);
}
}

#endif // FIND_CLOSEST_FACET_H
