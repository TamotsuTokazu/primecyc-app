#ifndef _RLWE_H_
#define _RLWE_H_

#include "binfhecontext.h"
#include <vector>

template <typename Poly>
struct SchemeImpl {
    using Vector = typename Poly::Vector;
    using Integer = typename Poly::Integer;
    using ILParams = typename Poly::Params;

    using RLWECiphertext = std::vector<Poly>;
    using RLWEKey = std::vector<Poly>;
    using RLWESwitchingKey = std::vector<std::vector<RLWECiphertext>>;
    using RGSWCiphertext = std::pair<std::vector<RLWECiphertext>, std::vector<RLWECiphertext>>;

    struct Params {
        std::shared_ptr<ILParams> poly;
        uint32_t p;
        Integer Q;
        uint32_t Bks;
    };

    Params params;

    std::vector<RLWESwitchingKey> ksk_galois;
    std::vector<RGSWCiphertext> bk;

    std::vector<uint32_t> iso_indices;

    Poly skp;

    SchemeImpl(Params p, Vector sk = {});

    Poly GaloisConjugate(const Poly &x, const uint32_t &a);

    template <typename T>
    std::vector<T> GaloisConjugate(const std::vector<T> &x, const uint32_t &a);

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain);
    Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain);

    void ModSwitch(Poly &x, const Integer &q);
    void ModSwitch(RLWECiphertext &ct, const Integer &q);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);

    RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);

    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);

    RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWECiphertext Process(Vector a, Integer b, Integer q_plain);
};

#include "rlwe-impl.h"

#endif