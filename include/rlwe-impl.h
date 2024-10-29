#ifndef _RLWE_IMPL_H_
#define _RLWE_IMPL_H_

#include "params.h"
#include "rlwe.h"

Scheme::Scheme(Params p, Vector sk) : params(p), ksk_galois(params.p), bk(sk.GetLength()), skp(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT) {
    sk.SetModulus(params.p);
    sk.ModEq(params.p);
    skp.SetFormat(EVALUATION);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(params.p - 2))
    for (uint32_t i = 2; i < params.p; i++) {
        auto skpi = GaloisConjugate(skp, i);
        auto ksk = KeySwitchGen({skpi}, {skp});
        ksk_galois[i] = ksk;
    }
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(sk.GetLength()))
    for (uint32_t i = 0; i < sk.GetLength(); i++) {
        Poly m(params.poly, COEFFICIENT, true);
        auto t = sk[i].ConvertToInt();
        if (t < params.p - 1) {
            m[sk[i].ConvertToInt()] = 1;
        } else {
            for (uint32_t j = 0; j < params.p - 1; j++) {
                m[j] = Integer(params.poly->GetModulus() - 1);
            }
        }
        m.SetFormat(EVALUATION);
        bk[i] = RGSWEncrypt(m, {skp});
    }
}

Poly Scheme::GaloisConjugate(const Poly &x, const uint32_t &a) {
    auto z = x;
    uint32_t n = x.GetLength() + 1;
    z.SetFormat(EVALUATION);
    Poly y(params.poly, EVALUATION, true);
    for (uint32_t i = 1; i < n; ++i) {
        y[i - 1] = z[(i * a) % n - 1];
    }
    return y;
}

template <typename T>
std::vector<T> Scheme::GaloisConjugate(const std::vector<T> &x, const uint32_t &a) {
    std::vector<T> y;
    for (auto &xi : x) {
        y.push_back(GaloisConjugate(xi, a));
    }
    return y;
}

RLWECiphertext Scheme::RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain) {
    uint32_t k = sk.size();
    lbcrypto::DiscreteUniformGeneratorImpl<Vector> dug;
    RLWECiphertext ct;

    Poly result(params.poly, EVALUATION, true);
    for (uint32_t i = 0; i < k; ++i) {
        Poly a(dug, params.poly, EVALUATION);
        result += a * sk[i];
        ct.push_back(std::move(a));
    }
    Poly e(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT);
    e.SetFormat(EVALUATION);
    result += e + m * (p::Q / q_plain);
    ct.push_back(std::move(result));
    return ct;
}

Poly Scheme::RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain) {
    uint32_t k = sk.size();
    Poly result = ct[k];
    for (uint32_t i = 0; i < k; ++i) {
        result -= ct[i] * sk[i];
    }
    result.SetFormat(COEFFICIENT);
    ModSwitch(result, q_plain);
    return result;
}

void Scheme::ModSwitch(Poly &x, const Integer &q) {
    uint32_t n = x.GetLength();
    Integer Q = x.GetModulus();
    Integer halfQ = Q / 2;
    for (uint32_t i = 0; i < n; ++i) {
        x[i] = (x[i] * q + halfQ) / Q;
        if (x[i] >= q) {
            x[i] -= q;
        }
    }
}

void Scheme::ModSwitch(RLWECiphertext &ct, const Integer &q) {
    uint32_t k = ct.size();
    for (uint32_t i = 0; i < k; ++i) {
        ModSwitch(ct[i], q);
    }
}

RLWESwitchingKey Scheme::KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN) {
    RLWESwitchingKey result;
    for (Integer t = 1; t <= p::Q; t *= p::Bks) {
        std::vector<RLWECiphertext> vec;
        for (auto &si : sk) {
            vec.push_back(RLWEEncrypt(t * si, skN, p::Q));
        }
        result.push_back(std::move(vec));
    }
    return result;
}

RLWECiphertext Scheme::KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K) {
    uint32_t k = ct.size() - 1;
    uint32_t kN = K[0][0].size() - 1;
    uint32_t l = K.size();
    RLWECiphertext result(kN + 1);
    for (uint32_t i = 0; i <= kN; ++i) {
        result[i] = Poly(params.poly, EVALUATION, true);
    }
    result[kN] += ct[k];
    for (uint32_t i = 0; i < k; ++i) {
        auto a = ct[i];
        a.SetFormat(COEFFICIENT);
        for (uint32_t j = 0; j < l; ++j) {
            Poly a0(params.poly, COEFFICIENT, true);
            for (uint32_t m = 0; m < a.GetLength(); ++m) {
                a0[m] = a[m] % p::Bks;
                a[m] /= p::Bks;
            }
            a0.SetFormat(EVALUATION);
            for (uint32_t m = 0; m <= kN; ++m) {
                result[m] -= a0 * K[j][i][m];
            }
        }
    }
    return result;
}

RGSWCiphertext Scheme::RGSWEncrypt(const Poly &m, const RLWEKey &sk) {
    if (sk.size() != 1) {
        throw std::runtime_error("RGSW encryption requires a single key");
    }
    Poly s = sk[0];
    std::vector<RLWECiphertext> vc, vcs;
    for (Integer t = 1; t <= p::Q; t *= p::Bks) {
        RLWECiphertext c = RLWEEncrypt(t * m, sk, p::Q);
        RLWECiphertext cs = RLWEEncrypt(t * m * s, sk, p::Q);
        vc.push_back(c);
        vcs.push_back(cs);
    }
    return std::make_pair(std::move(vcs), std::move(vc));
}

RLWECiphertext Scheme::ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW) {
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
            ai[j] = a[j] % p::Bks;
            bi[j] = b[j] % p::Bks;
            a[j] /= p::Bks;
            b[j] /= p::Bks;
        }
        ai.SetFormat(EVALUATION);
        bi.SetFormat(EVALUATION);
        ra += ai * ctGSW.first[i][0] + bi * ctGSW.second[i][0];
        rb += ai * ctGSW.first[i][1] + bi * ctGSW.second[i][1];
    }
    return {std::move(ra), std::move(rb)};
}

RLWECiphertext Scheme::Process(Vector a, Integer b) {
    a.SetModulus(params.p);
    a.ModEq(params.p);
    b.ModEq(params.p);
    Poly ca(params.poly, COEFFICIENT, true);
    Poly cb(params.poly, COEFFICIENT, true);
    if (b == params.p - 1) {
        for (uint32_t i = 0; i < params.p - 1; i++) {
            cb[i] = params.poly->GetModulus() - p::Q / p::t;
        }
    } else {
        cb[b.ConvertToInt()] = p::Q / p::t;
    }
    ca.SetFormat(EVALUATION);
    cb.SetFormat(EVALUATION);
    RLWECiphertext c{ca, cb};
    Integer t = 1;
    for (uint32_t i = 0; i < bk.size(); i++) {
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