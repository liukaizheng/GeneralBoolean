#ifndef ORIENTED_FACETS_H
#define ORIENTED_FACETS_H

#include <Eigen/Dense>

namespace bo {
template <typename DerivedF, typename DerivedE>
inline void orientedFacets(
        const Eigen::MatrixBase<DerivedF> &F,
        Eigen::PlainObjectBase<DerivedE> &E)
{
    E.resize(F.rows()*F.cols(),F.cols()-1);
      typedef typename DerivedE::Scalar EScalar;
      switch(F.cols())
      {
        case 4:
          E.block(0*F.rows(),0,F.rows(),1) = F.col(1).template cast<EScalar>();
          E.block(0*F.rows(),1,F.rows(),1) = F.col(3).template cast<EScalar>();
          E.block(0*F.rows(),2,F.rows(),1) = F.col(2).template cast<EScalar>();

          E.block(1*F.rows(),0,F.rows(),1) = F.col(0).template cast<EScalar>();
          E.block(1*F.rows(),1,F.rows(),1) = F.col(2).template cast<EScalar>();
          E.block(1*F.rows(),2,F.rows(),1) = F.col(3).template cast<EScalar>();

          E.block(2*F.rows(),0,F.rows(),1) = F.col(0).template cast<EScalar>();
          E.block(2*F.rows(),1,F.rows(),1) = F.col(3).template cast<EScalar>();
          E.block(2*F.rows(),2,F.rows(),1) = F.col(1).template cast<EScalar>();

          E.block(3*F.rows(),0,F.rows(),1) = F.col(0).template cast<EScalar>();
          E.block(3*F.rows(),1,F.rows(),1) = F.col(1).template cast<EScalar>();
          E.block(3*F.rows(),2,F.rows(),1) = F.col(2).template cast<EScalar>();
          return;
        case 3:
          E.block(0*F.rows(),0,F.rows(),1) = F.col(1).template cast<EScalar>();
          E.block(0*F.rows(),1,F.rows(),1) = F.col(2).template cast<EScalar>();
          E.block(1*F.rows(),0,F.rows(),1) = F.col(2).template cast<EScalar>();
          E.block(1*F.rows(),1,F.rows(),1) = F.col(0).template cast<EScalar>();
          E.block(2*F.rows(),0,F.rows(),1) = F.col(0).template cast<EScalar>();
          E.block(2*F.rows(),1,F.rows(),1) = F.col(1).template cast<EScalar>();
          return;
      }
}
}

#endif // ORIENTED_FACETS_H
