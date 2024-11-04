#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "binfhecontext.h"

using Poly = lbcrypto::NativePoly;
using Vector = typename Poly::Vector;
using Integer = typename Poly::Integer;
using ILParams = Poly::Params;

using RLWECiphertext = std::vector<Poly>;
using RLWEKey = std::vector<Poly>;
using RLWESwitchingKey = std::vector<std::vector<RLWECiphertext>>;
using RGSWCiphertext = std::pair<std::vector<RLWECiphertext>, std::vector<RLWECiphertext>>;

namespace p {
    const Integer t = 1 << 6;
    const uint32_t n = 600;
    const uint32_t p0 = 1153;
    const uint32_t p1 = 1297;
    const uint32_t pq = p0 * p1;
    const uint32_t Bks = 1 << 8;
    
    // const Integer Q = lbcrypto::FirstPrime<Integer>(56, p0 * p1 * (p0 - 1) * (p1 - 1));
    // const Integer rootOfUnity = lbcrypto::RootOfUnity<Integer>(p0 * p1 * (p0 - 1) * (p1 - 1), Q);
    // const Integer rootOfUnity0 = rootOfUnity.ModExp(p1 * (p1 - 1), Q);
    // const Integer rootOfUnity1 = rootOfUnity.ModExp(p0 * (p0 - 1), Q);
    const Integer Q = 36066736134770689;
    const Integer rootOfUnity = 4364918564594134;
    const Integer rootOfUnity0 = 28679241126083710;
    const Integer rootOfUnity1 = 9730598941305417;
    const Integer rootOfUnitypq = 6492451311899152;

    const Integer bigModulus = Q;
    const Integer bigRootOfUnity0 = rootOfUnity0;
    const Integer bigRootOfUnity1 = rootOfUnity1;
    const Integer bigRootOfUnitypq = rootOfUnitypq;

    const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, bigModulus, bigRootOfUnity0);
    const auto pp1 = std::make_shared<ILParams>(p1, Q, rootOfUnity1, bigModulus, bigRootOfUnity1);
    const auto ppq = std::make_shared<ILParams>(pq, Q, rootOfUnitypq, bigModulus, bigRootOfUnitypq);

}

struct Params {
    std::shared_ptr<ILParams> poly;
    uint32_t p;
};

#endif