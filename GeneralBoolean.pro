TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

EXTERNAL = $$PWD/external

BOOST = $$EXTERNAL/boost_1_68_0
CGAL = $$EXTERNAL/CGAL_4_13
GMP = $$EXTERNAL/gmp
Eigen = $$EXTERNAL/eigen

INCLUDEPATH += $$BOOST/include \
    $$CGAL/include \
    $$GMP/include \
    $$Eigen

QMAKE_CXXFLAGS += /bigobj
win32:CONFIG(release, debug|release) {
    LIBS += -L$$CGAL/lib -llibCGAL-vc140-mt-4.13 \
        -L$$GMP/lib -llibgmp-10 -llibmpfr-4
}
win32:CONFIG(debug, debug|release) {
    LIBS += -L$$CGAL/lib -llibCGAL-vc140-mt-gd-4.13 \
        -L$$GMP/lib -llibgmp-10 -llibmpfr-4
}

SOURCES += \
        main.cpp

HEADERS += \
    helpers/mesh_to_cagal_triangle_list.h \
    helpers/projected_cdt.h \
    helpers/insert_into_cdt.h \
    helpers/assign_scalar.h \
    helpers/unique_rows.h \
    helpers/list_to_matrix.h \
    helpers/read_obj.h \
    helpers/self_intersect_mesh.h \
    helpers/write_obj.h \
    helpers/sort_rows.h \
    helpers/sort_cols.h \
    helpers/oriented_facets.h \
    helpers/unique_edge_map.h \
    helpers/extract_manifold_patches.h \
    helpers/extract_cells.h \
    helpers/order_facets_around_edge.h \
    helpers/extract_triangle_adjacency.h \
    helpers/extract_components.h \
    helpers/vertex_face_adjacency.h \
    helpers/out_vertex.h \
    helpers/out_edge.h \
    helpers/out_facet.h \
    helpers/submesh_aabb_tree.h \
    helpers/find_closest_facet.h \
    helpers/unordered_pair.h \
    helpers/cell_adjacency.h \
    mesh_boolean_type.h \
    binary_winding_number_operations.h \
    mesh_boolean.h \
    helpers/remesh_self_intersections.h \
    helpers/propagate_winding_numbers.h \
    helpers/piecewise_constant_winding_number.h \
    helpers/resolve_duplicate_faces.h \
    helpers/remove_unreferenced_vertices.h \
    helpers/merge_matrices.h \
    helpers/remove_closed_points.h \
    helpers/extract_patch_adjacency.h \
    helpers/writepatches.h
