#ifndef _RLWE_H_
#define _RLWE_H_

#include "binfhecontext.h"
#include <vector>

template <typename Poly>
class SchemeImpl {
public:
    using Vector = typename Poly::Vector;
    using Integer = typename Poly::Integer;
    using ILParams = typename Poly::Params;

    using RLWEKey = std::vector<Poly>;
    using RLWECiphertext = std::vector<Poly>;
    using RLWEGadgetCiphertext = std::vector<RLWECiphertext>;
    using RLWESwitchingKey = std::vector<RLWEGadgetCiphertext>;
    using RGSWCiphertext = std::pair<RLWEGadgetCiphertext, RLWEGadgetCiphertext>;

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

    Poly skp;

    SchemeImpl(Params p, bool keygen = true);
    void GaloisKeyGen();

    Poly GaloisConjugate(const Poly &x, const usint &a);

    template <typename T>
    std::vector<T> GaloisConjugate(const std::vector<T> &x, const usint &a);

    template <typename T1, typename T2>
    std::pair<T1, T2> GaloisConjugate(const std::pair<T1, T2> &x, const usint &a);

    void ModSwitch(Poly &x, const Integer &q);
    void ModSwitch(RLWECiphertext &ct, const Integer &q);

    virtual RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain);
    RLWEGadgetCiphertext RLWEGadgetEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain);
    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);
    Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain);

    RLWECiphertext Mult(Poly a, RLWEGadgetCiphertext ct);
    RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);
    RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);
};

#include "rlwe-impl.h"

#endif