#ifndef REMESH_SELF_INTERSECTIONS_H
#define REMESH_SELF_INTERSECTIONS_H

#include "self_intersect_mesh.h"

namespace bo {
    template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedVV,
        typename DerivedFF,
        typename DerivedIF,
        typename DerivedJ,
        typename DerivedIM>
inline void remeshSelfIntersections(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        Eigen::PlainObjectBase<DerivedVV> &VV,
        Eigen::PlainObjectBase<DerivedFF> &FF,
        Eigen::PlainObjectBase<DerivedIF> &IF,
        Eigen::PlainObjectBase<DerivedJ> &J,
        Eigen::PlainObjectBase<DerivedIM> &IM)
{
    SelfIntersectMesh <
            CGAL::Epeck,
            DerivedV,
            DerivedF,
            DerivedVV,
            DerivedFF,
            DerivedIF,
            DerivedJ,
            DerivedIM> SIM(V, F, VV, FF, IF, J, IM);
}
}

#endif // REMESH_SELF_INTERSECTIONS_H
