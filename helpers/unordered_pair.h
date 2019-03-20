#ifndef UNORDERED_PAIR_H
#define UNORDERED_PAIR_H

#include <algorithm>

namespace bo {

template <typename Scalar>
struct unordered_pair
{
    unordered_pair(const Scalar &x, const Scalar &y) {data_[0] = x; data_[1] = y;}
    Scalar data_[2];
    friend inline bool operator ==(const unordered_pair<Scalar> &a, const unordered_pair<Scalar> &b)
    {
        return (a.data_[0] == b.data_[0] && a.data_[1] == b.data_[1]) || (a.data_[0] == b.data_[1] && a.data_[1] == b.data_[0]);
    }
    friend inline bool operator <(const unordered_pair<Scalar> &a, const unordered_pair<Scalar> &b)
    {
        Scalar ad[2] = {a.data_[0], a.data_[1]};
        Scalar bd[2] = {b.data_[0], b.data_[1]};
        if(ad[0] > ad[1]) std::swap(ad[0], ad[1]);
        if(bd[0] > bd[1]) std::swap(bd[0], bd[1]);
        return ad[0] < bd[0] || (ad[0] == bd[0] && ad[1] < bd[1]);
    }
};
}

#endif // UNORDERED_PAIR_H
