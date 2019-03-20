#ifndef UNIQUE_EDGE_MAP_H
#define UNIQUE_EDGE_MAP_H

#include "oriented_facets.h"
#include "sort_cols.h"
#include "unique_rows.h"

namespace bo {
template <
  typename DerivedF,
  typename DerivedE,
  typename DeriveduE,
  typename DerivedEMAP,
  typename uE2EType>
inline void uniqueEdgeMap(
  const Eigen::MatrixBase<DerivedF> &F,
  Eigen::PlainObjectBase<DerivedE> &E,
  Eigen::PlainObjectBase<DeriveduE> &uE,
  Eigen::PlainObjectBase<DerivedEMAP> &EMAP,
  std::vector<std::vector<uE2EType> > &uE2E)
{
    orientedFacets(F, E);
    DerivedE sortE, tuE;
    DerivedEMAP e2u;
    sortCols(E,true, sortE);
    uniqueRows(sortE, tuE, EMAP, e2u);
    const size_t &ne = tuE.rows();
    uE.resize(ne, E.cols());
    for(size_t i = 0; i < ne; i++)
        uE.row(i) = E.row(e2u(i));
    uE2E.resize(ne);
    for(auto &u : uE2E) u.reserve(2);
    for(size_t i = 0; i < E.rows(); i++)
        uE2E[EMAP(i)].emplace_back(i);
}
}

#endif // UNIQUE_EDGE_MAP_H
