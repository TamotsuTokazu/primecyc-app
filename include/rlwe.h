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
        usint p;
        Integer Q;
        usint Bks;
        std::vector<Integer> g;

        Params(std::shared_ptr<ILParams> poly_, usint p_, Integer Q_, usint Bks_) : poly(poly_), p(p_), Q(Q_), Bks(Bks_), g() {
            Integer pw = 1;
            lbcrypto::BigInteger bQ(Q.ToString());
            for (lbcrypto::BigInteger t = 1; t <= bQ; t *= Bks) {
                g.push_back(pw);
                pw.ModMulEq(Bks, Q);
            }
        }
    };

    Params params;

    std::vector<RLWESwitchingKey> ksk_galois;
    std::vector<RGSWCiphertext> bk;

    Poly skp;

    SchemeImpl(Params p, Vector sk = {});

    Poly GaloisConjugate(const Poly &x, const usint &a);

    template <typename T>
    std::vector<T> GaloisConjugate(const std::vector<T> &x, const usint &a);

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain);
    Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain);

    void ModSwitch(Poly &x, const Integer &q);
    void ModSwitch(RLWECiphertext &ct, const Integer &q);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);

    RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);

    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);

    RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWECiphertext Process(Vector a, Integer b, Integer q_plain, Integer mult);
};

#include "rlwe-impl.h"

#endif