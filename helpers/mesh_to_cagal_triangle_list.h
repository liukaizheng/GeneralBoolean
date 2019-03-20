#ifndef MESH_TO_CAGAL_TRIANGLE_LIST_H
#define MESH_TO_CAGAL_TRIANGLE_LIST_H

#include <vector>
#include <Eigen/Dense>
#include <CGAL/Triangle_3.h>

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename Kernel>
inline void meshToCGALTriangleList(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        std::vector<CGAL::Triangle_3<Kernel>> &T)
{
    typedef CGAL::Triangle_3<Kernel> Triangle_3;
    typedef typename Kernel::Point_3 Point_3;
    assert(V.cols() == 3);
    Eigen::Matrix<
            typename Kernel::FT,
            DerivedV::RowsAtCompileTime,
            DerivedV::ColsAtCompileTime>
            KV(V.rows(), V.cols());
    for(size_t i = 0; i < V.rows(); i++) {
        for(size_t j = 0; j < V.cols(); j++) {
            KV(i,j) = V(i,j);
        }
    }
    for(size_t i = 0; i < F.rows(); i++) {
        T.push_back(Triangle_3(
                        Point_3(KV(F(i,0), 0), KV(F(i,0), 1), KV(F(i,0), 2)),
                        Point_3(KV(F(i,1), 0), KV(F(i,1), 1), KV(F(i,1), 2)),
                        Point_3(KV(F(i,2), 0), KV(F(i,2), 1), KV(F(i,2), 2))));
    }
}
}
#endif // MESH_TO_CAGAL_TRIANGLE_LIST_H
