#ifndef PROJECTED_CDT_H
#define PROJECTED_CDT_H

#include "insert_into_cdt.h"
#include <unordered_map>

namespace bo {
template <typename Kernel, typename Index>
inline void projectedCDT(
        const std::vector<CGAL::Object> &objects,
        const typename Kernel::Plane_3 &P,
        std::vector<typename Kernel::Point_3> &vertices,
        std::vector<std::vector<Index>> &faces)
{
    typedef CGAL::Triangulation_vertex_base_2<Kernel>  TVB_2;
    typedef CGAL::Constrained_triangulation_face_base_2<Kernel> CTFB_2;
    typedef CGAL::Triangulation_data_structure_2<TVB_2,CTFB_2> TDS_2;
    typedef CGAL::Exact_intersections_tag Itag;
    typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS_2,Itag> CDT_2;
    typedef CGAL::Constrained_triangulation_plus_2<CDT_2> CDT_plus_2;
    CDT_plus_2 cdt;
    for(const auto &obj : objects) insert_into_cdt(obj, P, cdt);
    std::unordered_map<typename CDT_plus_2::Vertex_handle, Index> v2i;
    Index count = 0;
    for(auto it = cdt.finite_vertices_begin(); it != cdt.finite_vertices_end(); it++) {
        vertices.emplace_back(P.to_3d(it->point()));
        v2i[it] = count++;
    }
    for(auto it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); it++) {
        faces.push_back({v2i[it->vertex(0)], v2i[it->vertex(1)], v2i[it->vertex(2)]});
    }
}
}
#endif // PROJECTED_CDT_H
