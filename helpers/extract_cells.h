#ifndef EXTRACT_CELLS_H
#define EXTRACT_CELLS_H

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <algorithm>
#include "order_facets_around_edge.h"
#include "extract_triangle_adjacency.h"
#include "extract_components.h"
#include "out_facet.h"
#include "submesh_aabb_tree.h"
#include "find_closest_facet.h"
#include "vertex_face_adjacency.h"

namespace bo {
template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedP,
        typename DeriveduE,
        typename uE2EType,
        typename DerivedC>
inline typename DerivedF::Scalar extractCellsSingleComponent(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedP> &P,
        const Eigen::MatrixBase<DeriveduE> &uE,
        const std::vector<std::vector<uE2EType>> &uE2E,
        Eigen::PlainObjectBase<DerivedC> &cell)
{
    typedef typename DerivedF::Scalar Index;
    const Index &num_faces = F.rows();

    const auto edge_to_face = [&](const Index &ei) {return ei % num_faces;};

    const auto is_consistent =
       [&F](const size_t fid, const size_t s, const size_t d) -> bool {
       if ((size_t)F(fid, 0) == s && (size_t)F(fid, 1) == d) return false;
       if ((size_t)F(fid, 1) == s && (size_t)F(fid, 2) == d) return false;
       if ((size_t)F(fid, 2) == s && (size_t)F(fid, 0) == d) return false;

       if ((size_t)F(fid, 0) == d && (size_t)F(fid, 1) == s) return true;
       if ((size_t)F(fid, 1) == d && (size_t)F(fid, 2) == s) return true;
       if ((size_t)F(fid, 2) == d && (size_t)F(fid, 0) == s) return true;
       throw std::runtime_error("invalid faces");
       return false;
     };

    const Index &num_unique_edge = uE.rows();
    const Index &&num_patches = P.maxCoeff() + 1;

    std::vector<std::unordered_map<Index,Index>> patch_adj(num_patches);

    //这一步似乎非必要,下面那一步遍历uE instead of patch似乎就行
    for(Index i = 0; i < num_unique_edge; i++) {
        const auto &adj_faces = uE2E[i];
        const Index num_adj_faces = adj_faces.size();
        if(num_adj_faces > 2) {
            for(Index j = 0; j < num_adj_faces-1; j++) {
                const Index &patch_j = P(edge_to_face(adj_faces[j]));
                for(Index k = j+1; k < num_adj_faces; k++) {
                    const Index &patch_k = P(edge_to_face(adj_faces[k]));
                    if(patch_adj[patch_j].find(patch_k) == patch_adj[patch_j].end())
                        patch_adj[patch_j].emplace(patch_k, i);
                    if(patch_adj[patch_k].find(patch_j) == patch_adj[patch_k].end())
                        patch_adj[patch_k].emplace(patch_j, i);
                }
            }
        }
    }

    const Index INVALID = num_patches * 2;
    std::vector<std::unordered_set<Index>> equivalent_cells(num_patches * 2);
    std::vector<bool> processed(num_unique_edge, false);
    for(Index i = 0; i < num_patches; i++) {
        for(const auto &entry : patch_adj[i]) {
            const Index uei = entry.second;
            if(processed[uei]) continue;
            processed[uei] = true;

            const auto &adj_faces = uE2E[uei];
            const Index &num_adj_faces = adj_faces.size();
            assert(num_adj_faces > 2);

            const Index &s = uE(uei, 0);
            const Index &d = uE(uei, 1);
            std::vector<int> signed_adj_faces;
            for(const auto &ej : adj_faces) {
                const Index fid = edge_to_face(ej);
                bool cons = is_consistent(fid, s, d);
                signed_adj_faces.emplace_back(cons ? (fid + 1) : (-fid - 1));
            }

            Eigen::VectorXi order;
            orderFacetAroundEdge(V, F, s, d, signed_adj_faces, order);
            for(Index j = 0; j < num_adj_faces; j++) {
                const Index curr_idx = j;
                const Index next_idx = (j + 1) % num_adj_faces;
                const Index curr_patch_index = P(edge_to_face(adj_faces[order[curr_idx]]));
                const Index next_patch_index = P(edge_to_face(adj_faces[order[next_idx]]));
                const bool curr_cons = signed_adj_faces[order[curr_idx]] > 0;
                const bool next_cons = signed_adj_faces[order[next_idx]] > 0;
                const Index curr_cell = curr_patch_index * 2 + (curr_cons ? 0 : 1);
                const Index next_cell = next_patch_index * 2 + (next_cons ? 1 : 0);
                equivalent_cells[curr_cell].emplace(next_cell);
                equivalent_cells[next_cell].emplace(curr_cell);
            }
        }
    }

    Index count = 0;
    cell.resize(num_patches, 2);
    cell.setConstant(INVALID);
    auto extract_equivalent_cells = [&](const Index &i) {
        if(cell(i/2,i%2) != INVALID) return;
        std::queue<Index> Q;
        Q.push(i);
        cell(i/2, i%2) = count;
        while(!Q.empty()) {
             const Index &cur = Q.front();
             Q.pop();
             for(const auto &j : equivalent_cells[cur]) {
                 if(cell(j/2, j%2) != INVALID) continue;
                 Q.push(j);
                 cell(j/2, j%2) = count;
             }
        }
        count++;
    };
    for(Index i = 0; i < num_patches; i++) {
        extract_equivalent_cells(2 * i);
        extract_equivalent_cells(2 * i + 1);
    }
    assert((cell.array() != INVALID).all());
    return count;
}

template <
        typename DerivedV,
        typename DerivedF,
        typename DerivedP,
        typename DerivedE,
        typename DeriveduE,
        typename uE2EType,
        typename DerivedEMAP,
        typename DerivedC>
size_t extractCells(
        const Eigen::MatrixBase<DerivedV> &V,
        const Eigen::MatrixBase<DerivedF> &F,
        const Eigen::MatrixBase<DerivedP> &P,
        const Eigen::MatrixBase<DerivedE> &E,
        const Eigen::MatrixBase<DeriveduE> &uE,
        const std::vector<std::vector<uE2EType>> &uE2E,
        const Eigen::MatrixBase<DerivedEMAP> &EMAP,
        Eigen::PlainObjectBase<DerivedC> &cells)
{
    if(P.size() == 0) {
        assert(F.size() == 0);
        cells.resize(0, 2);
        return 0;
    }
    typedef typename DerivedF::Scalar Index;
    typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;

    const Eigen::Index &num_faces = F.rows();
    const size_t num_patches = P.maxCoeff() + 1;


    DerivedC raw_cells;
    const Index num_raw_cells = extractCellsSingleComponent(V, F, P, uE, uE2E, raw_cells);

    std::vector<std::vector<std::vector<size_t>>> TT;
    extractTriangleAdjcency(EMAP, uE2E, TT);

    DerivedP C;
    const size_t &num_components = extractComponents(TT, C);

    std::vector<std::vector<size_t>> components(num_components);
    for(Eigen::Index i = 0; i < num_faces; i++)
        components[C(i)].emplace_back(i);

    DerivedP outer_facets(num_components);
    DerivedP outer_cells(num_components);
    std::vector<DerivedP> Is(num_components);
    for(size_t i = 0; i < num_components; i++) {
        Is[i].resize(components[i].size());
        std::copy(components[i].begin(), components[i].end(), Is[i].data());
        bool flipped = false;
        outFacet(V, F, Is[i], outer_facets(i), flipped);
        outer_cells(i) = raw_cells(P(outer_facets(i)), flipped ? 1 : 0);
    }

    std::unordered_multimap<size_t, size_t> equivalent_cells;
    if(num_components > 1) {
        std::vector<
                CGAL::AABB_tree<
                    CGAL::AABB_traits<
                        Kernel,
                        CGAL::AABB_triangle_primitive<
                            Kernel,
                            std::vector<Kernel::Triangle_3 >::iterator >>>> trees(num_components);
        std::vector<std::vector<Kernel::Triangle_3>> triangle_lists(num_components);
        std::vector<std::vector<bool>> in_Is(num_components);

        for(size_t i = 0; i < num_components; i++)
            submeshAABBTree(V, F, Is[i], trees[i], triangle_lists[i], in_Is[i]);

        std::vector<std::vector<size_t>> VF, VFi;
        vertexFaceAdjacency(V.rows(), F, VF, VFi);

        DerivedV bbox_min(num_components, 3);
        DerivedV bbox_max(num_components, 3);
        bbox_min.setConstant(std::numeric_limits<double>::max());
        bbox_max.setConstant(std::numeric_limits<double>::min());

        for(Eigen::Index i = 0; i < num_faces; i++) {
            const auto &comp_id = C(i);
            const auto &f = F.row(i);
            for(size_t j = 0; j < 3; j++)
                for(size_t d = 0; d < 3; d++) {
                    bbox_min(comp_id, d) = std::min(bbox_min(comp_id,d), V(f(j),d));
                    bbox_max(comp_id, d) = std::max(bbox_max(comp_id,d), V(f(j),d));
                }
        }

        const auto bbox_intersects = [&bbox_max,&bbox_min](const size_t &ci, const size_t &cj)
        {
            return !(
                        bbox_max(ci,0) < bbox_min(cj,0) ||
                        bbox_max(ci,1) < bbox_min(cj,1) ||
                        bbox_max(ci,2) < bbox_min(cj,2) ||
                        bbox_max(cj,0) < bbox_min(ci,0) ||
                        bbox_max(cj,1) < bbox_min(ci,1) ||
                        bbox_max(cj,2) < bbox_min(ci,2));
        };

        const auto get_triangle_center = [&V,&F](const size_t &fid)
        {
            return ((V.row(F(fid,0))+V.row(F(fid,1))+V.row(F(fid,2)))/3.0).eval();
        };

        for(size_t i = 0; i < num_components; i++) {
            const auto &triangle_center = get_triangle_center(outer_facets(i));
            std::vector<std::pair<size_t, CGAL::Epeck::FT>> candidate_componets;
            for(size_t j = 0; j < num_components; j++) {
                if(i == j) continue;
                if(!bbox_intersects(i, j)) continue;
                size_t closet_face;
                bool is_flip;
                auto &&square_dist = findClosetFacet(V, F, Is[j], triangle_center, uE2E, EMAP, VF, VFi, trees[j], triangle_lists[j], in_Is[j], closet_face, is_flip);
                //ambient_cell 表示包含triangle_center的cell
                const size_t &ambient_cell = raw_cells(P(closet_face), is_flip ? 1 : 0);
                if(ambient_cell != outer_cells[j]) {
                    candidate_componets.emplace_back(ambient_cell, square_dist);
                }
            }
            if(candidate_componets.size()) {
                const auto &min_ele = std::min_element(candidate_componets.begin(), candidate_componets.end(), [](const std::pair<size_t, CGAL::Epeck::FT> &a, const std::pair<size_t, CGAL::Epeck::FT> &b) {return a < b;});
                equivalent_cells.emplace(outer_cells(i), min_ele->first);
                equivalent_cells.emplace(min_ele->first, outer_cells(i));
            }
        }
    }

    //所有等价的外部cell，接下来将会把它们都设为0
    std::unordered_set<size_t> outer_equivalent_cells;
    std::for_each(outer_cells.data(), outer_cells.data() + num_components, [&](const auto &oc) {
        const auto &it = equivalent_cells.find(oc);
        if(it == equivalent_cells.end())
            outer_equivalent_cells.emplace(oc);
    });

    //把外部cell都设为0,注意cells和raw_cells储存方式都必须为colmajor
    cells.resizeLike(raw_cells);
    const auto &INVALID = std::numeric_limits<typename DerivedC::Scalar>::max();
    cells.setConstant(INVALID);

    std::vector<size_t> new_cells(num_raw_cells, INVALID);
    for(const auto &idx : outer_equivalent_cells)
        new_cells[idx] = 0;

    typename DerivedC::Scalar count = 1;

    for(size_t i = 0; i < num_raw_cells; i++) {
        if(new_cells[i] != INVALID) continue;
        new_cells[i] = count;
        std::queue<size_t> Q;
        Q.emplace(i);
        while(!Q.empty()) {
            const auto &curr_cell = Q.front();
            Q.pop();
            const auto &e_pair = equivalent_cells.equal_range(curr_cell);
            std::for_each(e_pair.first, e_pair.second, [&](const auto &cell_pair) {
                const auto &equal_cell = cell_pair.second;
                if(new_cells[equal_cell] == INVALID) {
                    new_cells[equal_cell] = count;
                    Q.emplace(equal_cell);
                }
            });
        }
        count++;
    }

    for(size_t i = 0; i < num_patches; i++)
        for(size_t j = 0; j < 2; j++) {
            cells(i, j) = new_cells[raw_cells(i, j)];
        }

    assert((cells.array() != INVALID).all());
    return count;
}

}

#endif // EXTRACT_CELLS_H
