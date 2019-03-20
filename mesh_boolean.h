#ifndef MESH_BOOLEAN_H
#define MESH_BOOLEAN_H

#include "binary_winding_number_operations.h"
#include "helpers/remesh_self_intersections.h"
#include "helpers/unique_edge_map.h"
#include "helpers/extract_manifold_patches.h"
#include "helpers/extract_cells.h"
#include "helpers/propagate_winding_numbers.h"
#include "helpers/resolve_duplicate_faces.h"
#include "helpers/remove_unreferenced_vertices.h"
#include "helpers/merge_matrices.h"
#include "helpers/remove_boundary_patches.h"
#include "helpers/build_cell_mesh.h"
#include "helpers/extract_patch_adjacency.h"
#include "helpers/writepatches.h"
#include <functional>
#include <chrono>

#define print_msg(msg) \
    std::cout << #msg << "\t" << msg << "\n";

namespace bo {
template <
        typename DerivedVV,
        typename DerivedFF,
        typename Derivedsizes,
        typename DerivedVC,
        typename DerivedFC,
        typename DerivedJ>
bool meshBoolean(
        const Eigen::MatrixBase<DerivedVV> &VV,
        const Eigen::MatrixBase<DerivedFF> &FF,
        const Eigen::MatrixBase<Derivedsizes> &sizes,
        const std::function<bool (const Eigen::Matrix<int, Eigen::Dynamic, 1> &win_num)> &wind_num_op,
        Eigen::PlainObjectBase<DerivedVC> &VC,
        Eigen::PlainObjectBase<DerivedFC> &FC,
//        std::vector<DerivedVC> &Vs,
//        std::vector<DerivedFC> &Fs,
        Eigen::PlainObjectBase<DerivedJ> &J)
{
    typedef typename DerivedVC::Scalar Scalar;
    typedef CGAL::Epeck Kernel;
    typedef Kernel::FT ExactScalar;
    typedef Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic, 1> VectorXJ;
    typedef Eigen::Matrix<
            ExactScalar,
            Eigen::Dynamic,
            Eigen::Dynamic> MatrixXES;
    typedef Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXui;
    typedef Eigen::Matrix<size_t, Eigen::Dynamic, 1> VectorXui;

    MatrixXES V;
    DerivedFC F;
    VectorXJ CJ;
    auto start_time = std::chrono::high_resolution_clock::now();
    {
        Eigen::VectorXi I;
        Eigen::MatrixXi IF;
        remeshSelfIntersections(VV, FF, V, F, IF, CJ, I);
    }
    double spent_time = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start_time).count();
    std::cout << "resolved self-intersections spent " << spent_time / 1000<< " s\n";

    Eigen::MatrixXd dV;
    dV.resizeLike(V);
    for(Eigen::Index i = 0; i <dV.rows(); i++) {
        for(Eigen::Index j = 0; j < dV.cols(); j++)
            assignScalar(V(i, j), dV(i, j));
    }
    bo::writeOBJ("total.obj", dV, F);

    MatrixXui E, uE;
    VectorXui EMAP;
    std::vector<std::vector<size_t>> uE2E;
    uniqueEdgeMap(F, E, uE, EMAP, uE2E);

    VectorXui P;
    start_time = std::chrono::high_resolution_clock::now();
    extractManifoldPatches(F, EMAP, uE2E, P);
//    writePatches(V, F, P);

    MatrixXES rV;
    DerivedFC rF;
    VectorXJ rJ;
    VectorXui oP, nP;
    const auto num_old_patches = removeBoundaryPatches(V, F, P, CJ, uE2E, rV, rF, oP, rJ);
    uE2E.clear();
    uniqueEdgeMap(rF, E, uE, EMAP, uE2E);
    const auto num_new_patches = extractManifoldPatches(rF, EMAP, uE2E, nP);
    std::vector<unsigned> op2p(num_old_patches);
    std::vector<std::unordered_set<unsigned>> p2op(num_new_patches);
    if(num_old_patches == num_new_patches) {
        for(unsigned i = 0; i < num_old_patches; i++) {
            op2p[i] = i;
            p2op[i].emplace(i);
        }
    }
    else {
        for(Eigen::Index i = 0; i < rF.rows(); i++) {
            op2p[oP[i]] = nP[i];
            p2op[nP[i]].emplace(oP[i]);
        }
    }
    spent_time = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start_time).count();
    std::cout << "extract patches: " << spent_time / 1000<< " s\n";

    writePatches(rV,rF, oP);
//    std::vector<std::vector<Eigen::Index>> patch2face(num_patches);
//    for(Eigen::Index i = 0; i < rF.rows(); i++)
//        patch2face[nP(i)].emplace_back(i);


    std::vector<std::unordered_set<unsigned>> P2P;
    std::vector<std::unordered_map<unsigned, bool>> P2P_D;
    extractPatchAdjaceny(oP, E, uE2E, num_old_patches, P2P, P2P_D);


    MatrixXui per_patch_cells, PC;
    start_time = std::chrono::high_resolution_clock::now();
    const auto num_cells = extractCells(rV, rF, nP, E, uE, uE2E, EMAP, per_patch_cells);
    PC.resize(num_old_patches, 2);
    std::vector<std::unordered_set<unsigned>> cell2patch(num_cells);
    for(unsigned i = 0; i < num_new_patches; i++) {
        for(const auto& p: p2op[i]) {
            cell2patch[per_patch_cells(i, 0)].emplace(p);
            cell2patch[per_patch_cells(i, 1)].emplace(p);
            PC(p, 0) = per_patch_cells(i, 0);
            PC(p, 1) = per_patch_cells(i, 1);
        }
    }

    /*for(unsigned i = 0; i < num_cells; i++) {
        print_msg(i);
        for(const auto& p: cell2patch[i]) {
            print_msg(p);
        }
        std::cout << "\n";
    }*/

    Eigen::MatrixXi cW, W;
    VectorXui labels(rF.rows()), patch_labels(num_old_patches);
    const auto INVALID_PATCH_LABEL = std::numeric_limits<VectorXui::Scalar>::max();
    patch_labels.setConstant(INVALID_PATCH_LABEL);
    {
        std::vector<size_t> face_to_labels(FF.rows());
        size_t label_cout = 0, label_sum = sizes(0);
        for(size_t i = 0; i < FF.rows(); i++) {
            if(i >= label_sum) {
                label_cout++;
                label_sum += sizes(label_cout);
            }
            face_to_labels[i] = label_cout;
        }
        for(unsigned i = 0; i < rF.rows(); i++) {
            labels(i) = face_to_labels[rJ(i)];
            if(patch_labels(oP(i)) != INVALID_PATCH_LABEL)
                assert(patch_labels(oP(i)) == labels(i));
            else
                patch_labels(oP(i)) = labels(i);
        }

    }

    std::vector<std::unordered_map<unsigned, std::unordered_map<unsigned, bool>>> c2l2pd(num_cells); //cell -> label -> <patch,direction>
    for(unsigned i = 0; i < num_cells; i++) {
        std::unordered_map<unsigned, unsigned> L2P; //labels to patches
        std::unordered_set<unsigned> useless_labels;
        for(const auto& p: cell2patch[i]) {
                const auto& lab = patch_labels(p);
                const auto uit = useless_labels.find(lab);
                if(uit != useless_labels.end()) continue;

                const auto it = L2P.find(lab);
                if(it == L2P.end()) {
                    if(PC(p, 1) == i)
                        L2P.emplace(lab, p);
                    else
                        useless_labels.emplace(lab);
                }
        }

        print_msg(i);
        for(const auto& l2p: L2P) {
            const auto l= l2p.first;
            const auto bp = l2p.second; // patch for beginning
            std::unordered_map<unsigned, bool> non_flips(true);
            std::queue<unsigned> Q;
            Q.emplace(bp);
            while(!Q.empty()) {
                const auto p = Q.front();
                Q.pop();
                for(const auto& np: P2P[p]) {
                    const auto it = cell2patch[i].find(np);
                    if(it != cell2patch[i].end()) {
                        const auto lit = non_flips.find(np);
                        if(lit != non_flips.end()) continue;
                        Q.emplace(np);
                        bool non_flip = !(non_flips[p] ^ P2P_D[p][np]);
                        non_flips.emplace(np, non_flip);
                    }
                }
            }
            print_msg(l);
            for(const auto& it: non_flips) {
                std::cout << "p: " << it.first << "\t flip: " << it.second << "\n";
            }
            c2l2pd[i].emplace(l, std::move(non_flips));
        }
    }
    progateWindingNumbers(c2l2pd, rF, oP, patch_labels, PC, num_cells, cW, W);

    std::cout << cW << "\n";
    const auto num_faces = rF.rows();
    const size_t &num_inputs = sizes.size();
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> Wr(num_faces, 2);
    for(size_t i  = 0; i < num_faces; i++) {
        Eigen::VectorXi w_out(num_inputs), w_in(num_inputs);
        for(size_t j = 0; j < num_inputs; j++) {
            w_out(j) = W(i, 2*j);
            w_in(j) = W(i, 2*j+1);
        }
        Wr(i, 0) = wind_num_op(w_out);
        Wr(i, 1) = wind_num_op(w_in);
    }

    auto keep = [] (const bool &w_out, const bool &w_in) -> short {
        if(w_out && !w_in) return -1;
        else if(!w_out && w_in) return 1;
        else return 0;
    };

    auto index_to_signed_index = [] (const size_t &idx, bool ori) ->int {
        return int(idx + 1) * (ori ? 1 : -1);
    };

    std::vector<int> selected;
    for(size_t i = 0; i < num_faces; i++) {
        auto should_keep = keep(Wr(i,0), Wr(i,1));
        if(should_keep > 0)
            selected.emplace_back(index_to_signed_index(i, true));
        else if(should_keep < 0)
            selected.emplace_back(index_to_signed_index(i, false));
    }

    const size_t &num_selected = selected.size();
    DerivedFC keep_faces(num_selected, 3);
    DerivedJ keep_face_indices(num_selected);
    for(size_t i = 0; i < num_selected; i++) {
        const size_t idx= abs(selected[i]) - 1;
        if(selected[i] > 0)
            keep_faces.row(i) = rF.row(idx);
        else
            keep_faces.row(i) = rF.row(idx).reverse();

        keep_face_indices(i) = CJ(idx);
    }

    DerivedFC G;
    DerivedJ JJ;
    resolveDuplicateFaces(keep_faces, G, JJ);
    J.resizeLike(JJ);
    for(size_t i = 0; i < JJ.size(); i++)
        J(i) = keep_face_indices(JJ(i));

    Eigen::Matrix<typename DerivedVC::Scalar, Eigen::Dynamic, 3> Vs;
    Vs.resizeLike(rV);
    for(size_t i = 0; i < rV.rows(); i++)
        for(size_t j = 0; j < rV.cols(); j++)
            assignScalar(rV(i,j), Vs(i,j));

    removeUnreferencedVertices(Vs, G, VC, FC);

    return true;
}

template <
        typename DerivedV1,
        typename DerivedF1,
        typename DerivedV2,
        typename DerivedF2,
        typename DerivedVV,
        typename DerivedFF,
        typename DerivedJ>
bool meshBoolean(
        const Eigen::MatrixBase<DerivedV1> &V1,
        const Eigen::MatrixBase<DerivedF1> &F1,
        const Eigen::MatrixBase<DerivedV2> &V2,
        const Eigen::MatrixBase<DerivedF2> &F2,
        const MeshBooleanType &type,
        Eigen::PlainObjectBase<DerivedVV> &VV,
        Eigen::PlainObjectBase<DerivedFF> &FF,
//        std::vector<DerivedVV> &Vs,
//        std::vector<DerivedFF> &Fs,
        Eigen::PlainObjectBase<DerivedJ> &J)
{
    DerivedV1 mV;
    DerivedF1 mF;
    mergeMatrices(V1, F1, V2, F2, mV, mF);
    std::function<bool (const Eigen::Matrix<int, Eigen::Dynamic, 1> &)> win_num_op;
    meshBooleanTypeToFunction(type, win_num_op);
    DerivedJ sizes(2);
    sizes(0) = F1.rows();
    sizes(1) = F2.rows();
    return meshBoolean(mV, mF, sizes, win_num_op, VV, FF, J);
}
}


#endif // MESH_BOOLEAN_H
