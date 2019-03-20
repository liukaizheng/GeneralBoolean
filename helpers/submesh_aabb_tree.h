#ifndef SUBMESH_AABB_TREE_H
#define SUBMESH_AABB_TREE_H

#include <vector>
#include <Eigen/Dense>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedI,
        typename Kernel>
inline void submeshAABBTree(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedI> &I,
        CGAL::AABB_tree<
            CGAL::AABB_traits<
                Kernel,
                CGAL::AABB_triangle_primitive<
                    Kernel,
                    typename std::vector<typename Kernel::Triangle_3>::iterator>>> &tree,
        std::vector<typename Kernel::Triangle_3> &triangles,
        std::vector<bool> &in_I
        )
{
    in_I.resize(F.rows(), false);
    const Eigen::Index &num_faces = I.rows();
    for(Eigen::Index i = 0; i < num_faces; i++) {
        const auto &f = F.row(I(i));
        in_I[I(i)] = true;
        triangles.emplace_back(
                    typename Kernel::Point_3(V(f(0),0), V(f(0),1), V(f(0),2)),
                    typename Kernel::Point_3(V(f(1),0), V(f(1),1), V(f(1),2)),
                    typename Kernel::Point_3(V(f(2),0), V(f(2),1), V(f(2),2)));

        if(triangles.back().is_degenerate())
            assert(false && "degenerate triangles detected");
    }

    tree.insert(triangles.begin(), triangles.end());
    tree.accelerate_distance_queries();
}
}

#endif // SUBMESH_AABB_TREE_H
