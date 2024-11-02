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
    const uint32_t n = 2;
    // const uint32_t p0 = 1439;
    // const uint32_t p1 = 1447;
    const uint32_t p0 = 5;
    const uint32_t p1 = 7;
    const uint32_t pq = p0 * p1;
    const uint32_t Bks = 1 << 6;

    // const Integer Q("72057594241389101");
    // const Integer rootOfUnity0("3242320199475");
    // const Integer rootOfUnity1("1934345196961");
    // const Integer bigModulus("42535295865117307932921825928971042817");
    // const Integer bigRootOfUnity("34953128631915424799509995679592822");
    
    const Integer Q = lbcrypto::LastPrime<Integer>(56, p0 * p1 * (p0 - 1) * (p1 - 1) * t.ConvertToInt());
    // const Integer rootOfUnity0 = lbcrypto::RootOfUnity<Integer>(2 * p0, Q);
    // const Integer rootOfUnity1 = lbcrypto::RootOfUnity<Integer>(2 * p1, Q);
    const Integer rootOfUnity = lbcrypto::RootOfUnity<Integer>(p0 * p1 * (p0 - 1) * (p1 - 1), Q);
    const Integer rootOfUnity0 = rootOfUnity.ModExp(p1 * (p1 - 1), Q);
    const Integer rootOfUnity1 = rootOfUnity.ModExp(p0 * (p0 - 1), Q);
    const Integer rootOfUnitypq = lbcrypto::RootOfUnity<Integer>((p0 - 1) * (p1 - 1) / 2, Q);


    const uint32_t nttSize = 1 << (int)std::ceil(std::log2(2 * p0 - 1));
    const uint32_t nttSizepq = 1 << (int)std::ceil(std::log2(2 * pq - 1));
    const Integer bigModulus = lbcrypto::FirstPrime<Integer>(1 + (int)(std::ceil(std::log2(2 * pq - 1) + 2 * std::log2(Q.ConvertToInt()))), nttSizepq);
    const Integer bigRootOfUnity0 = lbcrypto::RootOfUnity<Integer>(nttSize, bigModulus);
    const Integer bigRootOfUnity1 = lbcrypto::RootOfUnity<Integer>(nttSize, bigModulus);
    const Integer bigRootOfUnitypq = lbcrypto::RootOfUnity<Integer>(nttSizepq, bigModulus);
    // const Integer bigModulus = Q;
    // const Integer bigRootOfUnity0 = rootOfUnity0;
    // const Integer bigRootOfUnity1 = rootOfUnity1;
    // const Integer bigRootOfUnitypq = rootOfUnitypq;

    const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, bigModulus, bigRootOfUnity0);
    const auto pp1 = std::make_shared<ILParams>(p1, Q, rootOfUnity1, bigModulus, bigRootOfUnity1);
    const auto ppq = std::make_shared<ILParams>(pq, Q, rootOfUnitypq, bigModulus, bigRootOfUnitypq);

}

struct Params {
    std::shared_ptr<ILParams> poly;
    uint32_t p;
};

#endif