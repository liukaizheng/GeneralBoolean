#ifndef ASSIGN_SCALAR_H
#define ASSIGN_SCALAR_H
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
namespace bo {
inline void assignScalar(const CGAL::Epeck::FT &_cgal, double &d)
{
    const CGAL::Epeck::FT &cgal = _cgal.exact();
    const auto &interval = CGAL::to_interval(cgal);
    d = interval.first;
    do {
        const double &next = nextafter(d, interval.second);
        if(CGAL::abs(d-cgal) < CGAL::abs(next - cgal)) break;
        d = next;
    }while (d < interval.second);
}

inline void assignScalar(const CGAL::Epeck::FT &_cgal, CGAL::Epeck::FT &d)
{
    d = _cgal;
}
}
#endif // ASSIGN_SCALAR_H
