#ifndef SELFINTERSECTMESH_H
#define SELFINTERSECTMESH_H

#include "mesh_to_cagal_triangle_list.h"
#include "projected_cdt.h"
#include "assign_scalar.h"
#include "unique_rows.h"
#include <CGAL/intersections.h>
#include <CGAL/box_intersection_d.h>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>

namespace bo {
template<
        typename Kernel,
        typename DerivedV,
        typename DerivedF,
        typename DerivedVV,
        typename DerivedFF,
        typename DerivedIF,
        typename DerivedJ,
        typename DerivedIM>
class SelfIntersectMesh
{
public:
    typedef SelfIntersectMesh<
    Kernel,
    DerivedV,
    DerivedF,
    DerivedVV,
    DerivedFF,
    DerivedIF,
    DerivedJ,
    DerivedIM> Self;

    typedef typename Kernel::Triangle_3 Triangle_3;
    typedef typename Kernel::Segment_3 Segment_3;
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Plane_3 Plane_3;
    typedef std::vector<Triangle_3> Triangles;
    typedef typename Triangles::iterator TrianglesIterator;
    typedef typename Triangles::const_iterator TrianglesConstIterator;
    typedef typename DerivedF::Scalar Index;
    typedef
        CGAL::Box_intersection_d::Box_with_handle_d<double,3,TrianglesIterator>
        Box;
    typedef std::vector<std::pair<Index, CGAL::Object>> ObjectList;

    SelfIntersectMesh(
            const Eigen::MatrixBase<DerivedV> &V,
            const Eigen::MatrixBase<DerivedF> &F,
            Eigen::PlainObjectBase<DerivedVV> &VV,
            Eigen::PlainObjectBase<DerivedFF> &FF,
            Eigen::PlainObjectBase<DerivedIF> &IF,
            Eigen::PlainObjectBase<DerivedJ> &J,
            Eigen::PlainObjectBase<DerivedIM> &IM);

    static void remeshIntersections(
            const Eigen::MatrixBase<DerivedV> &V,
            const Eigen::MatrixBase<DerivedF> &F,
            const std::vector<Triangle_3> &T,
            const std::unordered_map<Index, ObjectList> &offending,
            Eigen::PlainObjectBase<DerivedVV> &VV,
            Eigen::PlainObjectBase<DerivedFF> &FF,
            Eigen::PlainObjectBase<DerivedJ> &J,
            Eigen::PlainObjectBase<DerivedIM> &IM)
    {
        /*for(const auto &it : offending) {
            std::cout << it.first << ": ";
            for(const auto & i : it.second) {
                std::cout << i.first << " ";
            }
            std::cout << "\n";
        }*/
        const size_t &num_faces = F.rows();
        const size_t &num_base_vertices = V.rows();
        assert(num_faces == T.size());
        std::unordered_map<Index, std::vector<Index>> intersecting_and_coplanar;
        std::vector<bool> is_offending(num_faces, false);
        for(const auto &itr : offending) {
            const Index &fi = itr.first;
            is_offending[fi] = true;
            const auto &P = T[fi].supporting_plane();
            for(const auto &jtr: itr.second) {
                const auto &fj = jtr.first;
                const auto &tj = T[fj];
                if(P.has_on(tj[0]) && P.has_on(tj[1]) && P.has_on(tj[2])) {
                    const auto &loc = intersecting_and_coplanar.find(fi);
                    if(loc == intersecting_and_coplanar.end())
                        intersecting_and_coplanar[fi] = {fj};
                    else
                        loc->second.emplace_back(fj);
                }
            }
        }

        std::vector<std::vector<Index>> resolved_faces;
        std::vector<Index> source_faces;
        for(size_t i = 0; i < num_faces; i++) {
            if(!is_offending[i] && !T[i].is_degenerate()) {
                resolved_faces.push_back({F(i,0), F(i,1), F(i,2)});
                source_faces.emplace_back(i);
            }
        }

        std::vector<bool> processed(num_faces, false);
        std::vector<std::pair<Plane_3, std::vector<Index>>> cdt_inputs;
        for(const auto &itr : offending) {
            const auto &fi = itr.first;
            if(processed[fi]) continue;
            processed[fi] = true;
            std::vector<Index> involved_faces;
            const auto &loc = intersecting_and_coplanar.find(fi);
            if(loc == intersecting_and_coplanar.end())
                involved_faces.emplace_back(fi);
            else {
                std::queue<Index> Q;
                Q.push(fi);
                involved_faces.emplace_back(fi);
                while(Q.size() > 0) {
                    const auto &index = Q.front();
                    Q.pop();
                    const auto &overlapping_faces = intersecting_and_coplanar.find(index);
                    for(const auto &fj: overlapping_faces->second) {
                        if(processed[fj]) continue;
                        processed[fj] = true;
                        involved_faces.emplace_back(fj);
                        Q.push(fj);
                    }
                }
            }
            cdt_inputs.emplace_back(T[fi].supporting_plane(), involved_faces);
        }
        auto delaunay_triangulation = [&] (
                const Plane_3 &P,
                const std::vector<Index> &involved_faces,
                std::vector<Point_3> &vertices,
                std::vector<std::vector<Index>> &faces) -> void {
            std::vector<CGAL::Object> objects;
//            std::unordered_set<Index> processed_pairs;
            for(const Index &fid : involved_faces) {
                const auto &itr = offending.find(fid);
                const auto &t = T[fid];
                objects.emplace_back(CGAL::make_object(t));
                if(itr == offending.end()) {
                    continue;
                }else {
                    for(const auto & jtr: itr->second) {
//                        const auto fj = jtr.first;
//                        Index ind;
//                        if(fj < fid) ind = fj * num_faces + fid;
//                        else ind = fid * num_faces + fj;
//                        const auto &ind_it = processed_pairs.find(ind);
//                        if(ind_it == processed_pairs.end()) {
//                            processed_pairs.emplace(ind);
                            const auto &obj = jtr.second;
                            objects.emplace_back(obj);
//                        }
                    }
                }
            }
            projectedCDT<Kernel, Index>(objects, P, vertices, faces);
        };

        const size_t & num_cdt = cdt_inputs.size();
        std::vector<std::vector<Point_3>> cdt_vertices(num_cdt);
        std::vector<std::vector<std::vector<Index>>> cdt_faces(num_cdt);
        for(size_t i = 0; i < num_cdt; i++) {
            auto &vertices = cdt_vertices[i];
            auto &faces = cdt_faces[i];
            const auto &P = cdt_inputs[i].first;
            const auto &involved_faces = cdt_inputs[i].second;
            delaunay_triangulation(P, involved_faces, vertices, faces);
        }

        std::vector<Point_3> new_vertices;
        auto appendPoint = [&](const Point_3 &p) -> Index {
            const Index &ind = num_base_vertices + new_vertices.size();
            new_vertices.emplace_back(p);
            return ind;
        };
        auto postTriangulationProcess = [&](
                const std::vector<Point_3> &vertices,
                const std::vector<std::vector<Index>> &faces,
                const std::vector<Index> &involved_faces) {
            assert(involved_faces.size() > 0);
            // for all faces of the cdt
            for (const auto& f : faces)
            {
              const Point_3& v0 = vertices[f[0]];
              const Point_3& v1 = vertices[f[1]];
              const Point_3& v2 = vertices[f[2]];
              Point_3 center(
                (v0[0] + v1[0] + v2[0]) / 3.0,
                (v0[1] + v1[1] + v2[1]) / 3.0,
                (v0[2] + v1[2] + v2[2]) / 3.0);
              if (involved_faces.size() == 1)
              {
                // If only there is only one involved face, all sub-triangles must
                // belong to it and have the correct orientation.
                const auto& ori_f = involved_faces[0];
                std::vector<Index> corners(3);
                corners[0] = appendPoint(v0);
                corners[1] = appendPoint(v1);
                corners[2] = appendPoint(v2);
                resolved_faces.emplace_back(corners);
                source_faces.push_back(ori_f);
              } else
              {
                for (const auto& ori_f : involved_faces)
                {
                  const auto& triangle = T[ori_f];
                  const Plane_3 P = triangle.supporting_plane();
                  if (triangle.has_on(center)) {
                    std::vector<Index> corners(3);
                    corners[0] = appendPoint(v0);
                    corners[1] = appendPoint(v1);
                    corners[2] = appendPoint(v2);
                    if (CGAL::orientation(
                          P.to_2d(v0), P.to_2d(v1), P.to_2d(v2))
                        == CGAL::RIGHT_TURN) {
                      std::swap(corners[0], corners[1]);
                    }
                    resolved_faces.emplace_back(corners);
                    source_faces.push_back(ori_f);
                  }
                }
              }
            }
        };
        for(size_t i = 0; i < num_cdt; i++) {
            const auto &vertices = cdt_vertices[i];
            const auto &faces = cdt_faces[i];
            const auto &involved_faces = cdt_inputs[i].second;
            postTriangulationProcess(vertices, faces, involved_faces);
        }

        const size_t num_out_vertices = new_vertices.size() + num_base_vertices;
        VV.resize(num_out_vertices, 3);
        for(size_t i = 0; i < num_base_vertices; i++)
            for(size_t j = 0; j < 3; j++)
                assignScalar(V(i,j), VV(i,j));

        for(size_t i = 0; i < new_vertices.size(); i++)
            for(size_t j = 0; j < 3; j++)
                assignScalar(new_vertices[i][j], VV(num_base_vertices+i, j));

        const size_t &num_out_faces = resolved_faces.size();
        FF.resize(num_out_faces, 3);
        for(size_t i = 0; i < resolved_faces.size(); i++)
            for(size_t j = 0; j < 3; j++)
                FF(i, j) = resolved_faces[i][j];

        assert(source_faces.size() == num_out_faces);
        J.resize(num_out_faces);
        std::copy(source_faces.begin(), source_faces.end(), J.data());

        DerivedVV unique_vv;
        DerivedJ unique_to_vv, vv_to_unique;
        uniqueRows(VV, unique_vv, vv_to_unique, unique_to_vv);

        VV = unique_vv;

        std::transform(FF.data(), FF.data() + FF.rows()*FF.cols(), FF.data(), [&](const Index &a) {
            return vv_to_unique(a);
        });
        IM.resize(VV.rows());
        for(size_t i = 0; i < VV.rows(); i++)
            IM(i) = i;
    }

private:
    void processIntersectingBoxs(const std::vector<Triangle_3> & T, const std::vector<std::pair<TrianglesIterator, TrianglesIterator>> &candidate_triangle_list)
    {
        std::for_each(candidate_triangle_list.begin(), candidate_triangle_list.end(),[&](const std::pair<TrianglesIterator, TrianglesIterator> &tp) {
           Index fa = tp.first - T.begin();
           Index fb = tp.second - T.begin();
           const Triangle_3 &A = T[fa];
           const Triangle_3 &B = T[fb];
           Index total_shared_vertices = 0;
           std::vector<std::pair<Index, Index>> shared;
           for(Index i = 0; i < 3; i++)
               for(Index j = 0; j < 3; j++) {
                   if(F_(fa, i) == F_(fb, j) || A.vertex(i) == B.vertex(j)) {
                       total_shared_vertices++;
                       shared.emplace_back(i, j);
                   }
               }
           if(total_shared_vertices == 3)
               return;
           else if(total_shared_vertices == 2)
               doubleSharedVertex(A, B, fa, fb, shared);
           else if(total_shared_vertices == 1)
               singleSharedVertex(A, fa, shared[0].first, B, fb, shared[0].second);
           else
               intersect(A, fa, B, fb);
        });
    }

    bool doubleSharedVertex(
            const Triangle_3 &A,
            const Triangle_3 &B,
            const Index &ia,
            const Index &ib,
            const std::vector<std::pair<Index, Index>> &shared)
    {

        if(
                A.supporting_plane() != B.supporting_plane() &&
                A.supporting_plane() != B.supporting_plane().opposite())
            return false;

        auto oppositePointInside = [&](const Triangle_3 &A, const Index &i, const Index &j, const Triangle_3 &B) -> bool {
            const Index &o = opposite_index_[i+j-1];
            return CGAL::do_intersect(A.vertex(o), B);
        };
        auto oppositeEdgeIntersect = [&](const Triangle_3 &A, const Index &i, const Triangle_3 &B, const Index &j) -> bool {
            Segment_3 a(A.vertex((i+1)%3), A.vertex((i+2)%3));
            Segment_3 b(B.vertex((j+1)%3), B.vertex((j+2)%3));
            return CGAL::do_intersect(a, b);
        };
        if(
                !oppositePointInside(A, shared[0].first, shared[1].first, B) &&
                !oppositePointInside(B, shared[0].second, shared[1].second, A) &&
                !oppositeEdgeIntersect(A, shared[0].first, B, shared[1].second) &&
                !oppositeEdgeIntersect(A, shared[1].first, B, shared[0].second))
            return false;
        countIntersection(ia, ib);
        CGAL::Object result = CGAL::intersection(A, B);
        if(!result.is_empty()) {
            if(CGAL::object_cast<Segment_3>(&result) || CGAL::object_cast<Point_3>(&result)) {
                assert(false &&
                       "Co-planar non-degenerate triangles should intersect over triangle\n");
                return false;
            }else {
                offending_[ia].emplace_back(ib, result);
                offending_[ib].emplace_back(ia, result);
                return true;
            }
        }
        else {
            assert(false &&
                   "CGAL::intersection should agree with predicate test\n");
            return false;
        }
        return true;
    }

    bool singleSharedVertex(
            const Triangle_3 &A,
            const Index &fa,
            const Index &va,
            const Triangle_3 &B,
            const Index &fb,
            const Index &vb)
    {

        if(singleSharedVertex(A, fa, va, B, fb))
            return true;
        return singleSharedVertex(B, fb, vb, A, fa);
    }

    bool singleSharedVertex(
            const Triangle_3 &A,
            const Index &fa,
            const Index &va,
            const Triangle_3 &B,
            const Index &fb)
    {

        Segment_3 sa(A.vertex((va+1)%3), A.vertex((va+2)%3));
        if(!CGAL::do_intersect(sa, B)) return false;

        CGAL::Object result = CGAL::intersection(sa, B);
        if(const Point_3 *p = CGAL::object_cast<Point_3>(&result)) {
            CGAL::Object seg = CGAL::make_object(Segment_3(A.vertex(va), *p));
            countIntersection(fa, fb);
            offending_[fa].emplace_back(fb, seg);
            offending_[fb].emplace_back(fa, seg);
            return true;
        } else if (CGAL::object_cast<Segment_3>(&result)) {
            bool ret = intersect(A, fa, B, fb);
            if(!ret) {
                assert(false &&
                       "intersect should agree with do_intersect\n");
                return false;
            }
            return true;
        }
        return true;
    }

    bool intersect(
            const Triangle_3 &A,
            const Index &fa,
            const Triangle_3 &B,
            const Index &fb)
    {

        if(!CGAL::do_intersect(A, B)) return false;
        countIntersection(fa, fb);
        CGAL::Object result = CGAL::intersection(A, B);
        offending_[fa].emplace_back(fb, result);
        offending_[fb].emplace_back(fa, result);
        return true;
    }

    void countIntersection(const Index &fa, const Index &fb)
    {
        count_++;
        markOffensive(fa);
        markOffensive(fb);
    }

    void markOffensive(const Index &f)
    {

        lIF_.emplace_back(f);
        if(!offending_.count(f))
            offending_[f] = {};
    }

private:
    const Eigen::MatrixBase<DerivedV> &V_;
    const Eigen::MatrixBase<DerivedF> &F_;
    const Index opposite_index_[3] ={2, 1, 0};
    std::vector<Index> lIF_;
    std::unordered_map<Index, ObjectList> offending_;
    Index count_ = 0;
};

template<typename Kernel, typename DerivedV, typename DerivedF, typename DerivedVV, typename DerivedFF, typename DerivedIF, typename DerivedJ, typename DerivedIM>
inline SelfIntersectMesh<
    Kernel,
    DerivedV,
    DerivedF,
    DerivedVV,
    DerivedFF,
    DerivedIF,
    DerivedJ,
    DerivedIM>::SelfIntersectMesh(const Eigen::MatrixBase<DerivedV> &V, const Eigen::MatrixBase<DerivedF> &F, Eigen::PlainObjectBase<DerivedVV> &VV, Eigen::PlainObjectBase<DerivedFF> &FF, Eigen::PlainObjectBase<DerivedIF> &IF, Eigen::PlainObjectBase<DerivedJ> &J, Eigen::PlainObjectBase<DerivedIM> &IM):
    V_(V), F_(F), lIF_()
{
    std::vector<Triangle_3> triangle_list;
    meshToCGALTriangleList(V, F, triangle_list);
    std::vector<Box> boxs;
    boxs.reserve(F.rows());
    for(auto it = triangle_list.begin(); it != triangle_list.end(); it++) {
        if(!it->is_degenerate())
            boxs.emplace_back(Box(it->bbox(), it));
    }
    /*typedef Eigen::Matrix<Box::NT, Eigen::Dynamic, 3> M;
    M lows(boxs.size(), 3);
    M highs(boxs.size(), 3);
    for(unsigned i = 0; i < boxs.size(); i++) {
        for(unsigned j = 0; j < 3; j++) {
            lows(i, j) = boxs[i].min_coord(j);
            highs(i, j) = boxs[i].max_coord(j);
        }
    }
    std::cout << "l:\n" << lows;
    std::cout << "\n h: \n" << highs << "\n";*/
    std::vector<std::pair<TrianglesIterator, TrianglesIterator>> candidate_triangle_pairs;
    CGAL::box_self_intersection_d(boxs.begin(), boxs.end(), [&](const Box &a, const Box &b) {
        candidate_triangle_pairs.emplace_back(std::make_pair(a.handle(), b.handle()));
    });
    processIntersectingBoxs(triangle_list, candidate_triangle_pairs);
    assert(lIF_.size() % 2 == 0);
    IF.resize(lIF_.size()/ 2, 2);
    for(size_t i = 0; i < lIF_.size() / 2; i++) {
        IF(i, 0) = lIF_[2*i + 0];
        IF(i, 1) = lIF_[2*i + 1];
    }
    remeshIntersections(V, F, triangle_list, offending_, VV, FF, J, IM);
}

}

#endif
