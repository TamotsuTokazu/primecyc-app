#ifndef _RLWE_IMPL_H_
#define _RLWE_IMPL_H_

#include "rlwe.h"



template <typename Poly>
SchemeImpl<Poly>::SchemeImpl(Params p, Vector sk) : params(p), ksk_galois(params.p), bk(sk.GetLength()) {
    if (sk.GetLength() != 0) {
        skp = Poly(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT);
        sk.SetModulus(params.p);
        sk.ModEq(params.p);
        skp.SetFormat(EVALUATION);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(params.p - 2))
        for (usint i = 2; i < params.p; i++) {
            auto skpi = GaloisConjugate(skp, i);
            auto ksk = KeySwitchGen({skpi}, {skp});
            ksk_galois[i] = ksk;
        }
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(sk.GetLength()))
        for (usint i = 0; i < sk.GetLength(); i++) {
            Poly m(params.poly, COEFFICIENT, true);
            auto t = sk[i].ConvertToInt();
            if (t < params.p - 1) {
                m[sk[i].ConvertToInt()] = 1;
            } else {
                for (usint j = 0; j < params.p - 1; j++) {
                    m[j] = Integer(params.poly->GetModulus() - 1);
                }
            }
            m.SetFormat(EVALUATION);
            bk[i] = RGSWEncrypt(m, {skp});
        }
    }
}

template <typename Poly>
Poly SchemeImpl<Poly>::GaloisConjugate(const Poly &x, const usint &a) {
    if (x.GetFormat() != EVALUATION) {
        throw std::runtime_error("GaloisConjugate requires an evaluation format polynomial");
    }
    usint n = x.GetLength() + 1;
    Poly y(params.poly, EVALUATION, true);
    for (usint i = 1; i < n; ++i) {
        y[i - 1] = x[i * a % n - 1];
    }
    return y;
}

template <typename Poly> template <typename T>
std::vector<T> SchemeImpl<Poly>::GaloisConjugate(const std::vector<T> &x, const usint &a) {
    std::vector<T> y;
    for (auto &xi : x) {
        y.push_back(GaloisConjugate(xi, a));
    }
    return y;
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain) {
    usint k = sk.size();
    lbcrypto::DiscreteUniformGeneratorImpl<Vector> dug;
    RLWECiphertext ct;

    Poly result(params.poly, EVALUATION, true);
    for (usint i = 0; i < k; ++i) {
        Poly a(dug, params.poly, EVALUATION);
        result += a * sk[i];
        ct.push_back(std::move(a));
    }
    Poly e(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT);
    e.SetFormat(EVALUATION);
    result += e + m * (params.Q / q_plain);
    ct.push_back(std::move(result));
    return ct;
}

template <typename Poly>
Poly SchemeImpl<Poly>::RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain) {
    usint k = sk.size();
    Poly result = ct[k];
    for (usint i = 0; i < k; ++i) {
        result -= ct[i] * sk[i];
    }
    result.SetFormat(COEFFICIENT);
    ModSwitch(result, q_plain);
    return result;
}

template <typename Poly>
void SchemeImpl<Poly>::ModSwitch(Poly &x, const Integer &q) {
    using lbcrypto::BigInteger;
    usint n = x.GetLength();
    BigInteger Q = x.GetModulus();
    BigInteger halfQ = Q / 2;
    BigInteger qq = q;
    for (usint i = 0; i < n; ++i) {
        BigInteger xi = x[i];
        xi = (xi * qq + halfQ) / Q;
        if (xi >= qq) {
            xi -= qq;
        }
        x[i] = xi.ToString();
    }
}

template <typename Poly>
void SchemeImpl<Poly>::ModSwitch(RLWECiphertext &ct, const Integer &q) {
    usint k = ct.size();
    for (usint i = 0; i < k; ++i) {
        ModSwitch(ct[i], q);
    }
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWESwitchingKey SchemeImpl<Poly>::KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN) {
    RLWESwitchingKey result;
    for (Integer t = 1; t <= params.Q; t *= params.Bks) {
        std::vector<RLWECiphertext> vec;
        for (auto &si : sk) {
            vec.push_back(RLWEEncrypt(t * si, skN, params.Q));
        }
        result.push_back(std::move(vec));
    }
    return result;
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K) {
    usint k = ct.size() - 1;
    usint kN = K[0][0].size() - 1;
    usint l = K.size();
    RLWECiphertext result(kN + 1);
    for (usint i = 0; i <= kN; ++i) {
        result[i] = Poly(params.poly, EVALUATION, true);
    }
    result[kN] += ct[k];
    for (usint i = 0; i < k; ++i) {
        auto a = ct[i];
        a.SetFormat(COEFFICIENT);
        for (usint j = 0; j < l; ++j) {
            Poly a0(params.poly, COEFFICIENT, true);
            for (usint m = 0; m < a.GetLength(); ++m) {
                a0[m] = a[m] % params.Bks;
                a[m] /= params.Bks;
            }
            a0.SetFormat(EVALUATION);
            for (usint m = 0; m <= kN; ++m) {
                result[m] -= a0 * K[j][i][m];
            }
        }
    }
    return result;
}

template <typename Poly>
typename SchemeImpl<Poly>::RGSWCiphertext SchemeImpl<Poly>::RGSWEncrypt(const Poly &m, const RLWEKey &sk) {
    if (sk.size() != 1) {
        throw std::runtime_error("RGSW encryption requires a single key");
    }
    Poly s = sk[0];
    std::vector<RLWECiphertext> vc, vcs;
    for (Integer t = 1; t <= params.Q; t *= params.Bks) {
        RLWECiphertext c = RLWEEncrypt(t * m, sk, params.Q);
        RLWECiphertext cs = RLWEEncrypt(t * m * s, sk, params.Q);
        vc.push_back(c);
        vcs.push_back(cs);
    }
    return std::make_pair(std::move(vcs), std::move(vc));
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW) {
    if (ct.size() != 2) {
        throw std::runtime_error("ExtMult requires a ciphertext with two components");
    }
    Poly ra(params.poly, EVALUATION, true);
    Poly rb(params.poly, EVALUATION, true);
    u_int32_t l = ctGSW.first.size();
    Poly a = ct[0].Negate(), b = ct[1];
    a.SetFormat(COEFFICIENT);
    b.SetFormat(COEFFICIENT);
    for (u_int32_t i = 0; i < l; ++i) {
        Poly ai(params.poly, COEFFICIENT, true);
        Poly bi(params.poly, COEFFICIENT, true);
        for (u_int32_t j = 0; j < a.GetLength(); ++j) {
            ai[j] = a[j] % params.Bks;
            bi[j] = b[j] % params.Bks;
            a[j] /= params.Bks;
            b[j] /= params.Bks;
        }
        ai.SetFormat(EVALUATION);
        bi.SetFormat(EVALUATION);
        ra += ai * ctGSW.first[i][0] + bi * ctGSW.second[i][0];
        rb += ai * ctGSW.first[i][1] + bi * ctGSW.second[i][1];
    }
    return {std::move(ra), std::move(rb)};
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::Process(Vector a, Integer b, Integer q_plain) {
    a.SetModulus(params.p);
    a.ModEq(params.p);
    b.ModEq(params.p);
    Poly ca(params.poly, COEFFICIENT, true);
    Poly cb(params.poly, COEFFICIENT, true);
    if (b == params.p - 1) {
        for (usint i = 0; i < params.p - 1; i++) {
            cb[i] = params.poly->GetModulus() - params.Q / q_plain;
        }
    } else {
        cb[b.ConvertToInt()] = params.Q / q_plain;
    }
    ca.SetFormat(EVALUATION);
    cb.SetFormat(EVALUATION);
    RLWECiphertext c{ca, cb};
    Integer t = 1;
    for (usint i = 0; i < bk.size(); i++) {
        a[i] = params.p - a[i];
        if (a[i] != params.p) {
            t.ModMulEq(a[i].ModInverse(params.p), params.p);
            if (t != 1) {
                c = GaloisConjugate(c, t.ConvertToInt());
                c = KeySwitch(c, ksk_galois[t.ConvertToInt()]);
            }
            c = ExtMult(c, bk[i]);
            t = a[i];
        }
    }
    if (t != 1) {
        c = GaloisConjugate(c, t.ConvertToInt());
        c = KeySwitch(c, ksk_galois[t.ConvertToInt()]);
    }
    return c;
}

#endif