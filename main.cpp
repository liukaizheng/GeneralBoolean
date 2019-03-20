#include <iostream>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include "helpers/read_obj.h"
#include "helpers/write_obj.h"
#include "helpers/self_intersect_mesh.h"
#include "helpers/unique_edge_map.h"
#include "helpers/extract_manifold_patches.h"
#include "helpers/extract_cells.h"
#include "helpers/extrude_meshes.h"
#include "mesh_boolean.h"
#include "binary_winding_number_operations.h"
#include "helpers/resize_meshes.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
int main()
{
    typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXui;
    typedef Eigen::Matrix<size_t, Eigen::Dynamic, 1> VectorXui;
    Eigen::MatrixXd V1, V2, VV;
    MatrixXui F1, F2, FF;

    bo::readOBJ("./bp.obj", V1, F1);
    bo::readOBJ("./tp.obj", V2, F2);

    VectorXui J;
    bo::meshBoolean(V1, F1, V2, F2, bo::MeshBooleanType::MESH_BOOLEAN_INTERSECT, VV, FF, J);
    bo::writeOBJ("123.obj", VV, FF);
    return EXIT_SUCCESS;
}
