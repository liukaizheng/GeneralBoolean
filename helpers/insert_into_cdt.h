#ifndef INSERT_INTO_CDT_H
#define INSERT_INTO_CDT_H

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <vector>

namespace bo {
template <typename Kernel>
inline void insert_into_cdt(
        const CGAL::Object &object,
        const typename Kernel::Plane_3 &P,
        CGAL::Constrained_triangulation_plus_2<
            CGAL::Constrained_Delaunay_triangulation_2<
                Kernel,
                CGAL::Triangulation_data_structure_2<
                    CGAL::Triangulation_vertex_base_2<Kernel>,
                    CGAL::Constrained_triangulation_face_base_2<Kernel>
                >,
                CGAL::Exact_intersections_tag
            >
        > &cdt)
{
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Segment_3 Segment_3;
    typedef typename Kernel::Triangle_3 Triangle_3;
    if(const Segment_3 *seg = CGAL::object_cast<Segment_3>(&object))
        cdt.insert_constraint(P.to_2d(seg->vertex(0)), P.to_2d(seg->vertex(1)));
    else if(const Point_3 *point = CGAL::object_cast<Point_3>(&object))
        cdt.insert(P.to_2d(*point));
    else if(const Triangle_3 *triangle = CGAL::object_cast<Triangle_3>(&object)) {
        cdt.insert_constraint(P.to_2d(triangle->vertex(0)), P.to_2d(triangle->vertex(1)));
        cdt.insert_constraint(P.to_2d(triangle->vertex(1)), P.to_2d(triangle->vertex(2)));
        cdt.insert_constraint(P.to_2d(triangle->vertex(2)), P.to_2d(triangle->vertex(0)));
    }else if(const std::vector<Point_3> *poly = CGAL::object_cast<std::vector<Point_3>>(&object)) {
        const size_t & m = poly->size();
        assert(m > 2);
        for(size_t p = 0; p < m; p++) {
            const size_t np = (p+1) % m;
            cdt.insert_constraint(P.to_2d(poly->at(p)), P.to_2d(poly->at(np)));
        }
    }else
        assert(false && "error on inserting into constrained delaunay triangulation\n");
}
}


#endif // INSERT_INTO_CDT_H
