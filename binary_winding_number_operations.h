#ifndef BINARY_WINDING_NUMBER_OPERATIONS_H
#define BINARY_WINDING_NUMBER_OPERATIONS_H

#include "mesh_boolean_type.h"
#include <Eigen/Dense>

namespace bo {
template <MeshBooleanType Op>
struct BinaryWindingNumberOperations
{
    template<typename DerivedW>
    bool operator() (
            const Eigen::MatrixBase<DerivedW> &win_nums) const {
        assert(false && "not implemented");
        return false;
    }
};

template <>
struct BinaryWindingNumberOperations<MESH_BOOLEAN_UNION>
{
    template<typename DerivedW>
    bool operator () (
            const Eigen::MatrixBase<DerivedW> &win_nums) const {
        for(Eigen::Index i = 0; i < win_nums.size(); i++)
            if(win_nums(i) > 0) return true;
        return false;
    }
};

template <>
struct BinaryWindingNumberOperations<MESH_BOOLEAN_INTERSECT>
{
    template<typename DerivedW>
    bool operator () (
            const Eigen::MatrixBase<DerivedW> &win_nums) const {
        for(Eigen::Index i = 0; i < win_nums.size(); i++)
            if(win_nums(i) <= 0) return false;
        return true;
    }
};
template <>
struct BinaryWindingNumberOperations<MESH_BOOLEAN_MINUS>
{
    template<typename DerivedW>
    bool operator () (
            const Eigen::MatrixBase<DerivedW> &win_nums) const {
        assert(win_nums.size() > 1);
        bool union_rest = false;
        for(Eigen::Index i = 1; i < win_nums.size(); i++) {
            union_rest = union_rest || win_nums(i) > 0;
            if(union_rest) break;
        }
        return win_nums(0) > 0 && !union_rest;
    }
};

typedef BinaryWindingNumberOperations<MESH_BOOLEAN_UNION> BinaryUnion;
typedef BinaryWindingNumberOperations<MESH_BOOLEAN_INTERSECT> BinaryIntersect;
typedef BinaryWindingNumberOperations<MESH_BOOLEAN_MINUS> BinaryMinus;

inline void meshBooleanTypeToFunction(
        const MeshBooleanType &type,
        std::function<bool (const Eigen::Matrix<int, Eigen::Dynamic, 1> &)> &win_num_op)
{
    switch (type) {
    case MESH_BOOLEAN_UNION:
        win_num_op = BinaryUnion();
        break;
    case MESH_BOOLEAN_INTERSECT:
        win_num_op = BinaryIntersect();
        break;
    case MESH_BOOLEAN_MINUS:
        win_num_op = BinaryMinus();
        break;
    default:
        assert(false && "unsuported boolena type");
        break;
    }
}
}

#endif // BINARY_WINDING_NUMBER_OPERATIONS_H
