#ifndef _RLWE_IMPL_H_
#define _RLWE_IMPL_H_

template <typename Poly>
SchemeImpl<Poly>::SchemeImpl(Params p, bool keygen) : params(p) {
    if (keygen) {
        skp = Poly(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT);
    }
}

template <typename Poly>
void SchemeImpl<Poly>::GaloisKeyGen() {
    ksk_galois.resize(params.p);
    skp.SetFormat(EVALUATION);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(params.p - 2))
    for (usint i = 2; i < params.p; i++) {
        auto skpi = GaloisConjugate(skp, i);
        auto ksk = KeySwitchGen({skpi}, {skp});
        ksk_galois[i] = ksk;
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
        a *= 0;
        result += a * sk[i];
        ct.push_back(std::move(a));
    }
    Poly e(lbcrypto::DiscreteGaussianGeneratorImpl<Vector>(), params.poly, COEFFICIENT);
    e.SetFormat(EVALUATION);
    e *= 0;
    result += e + m * (params.Q / q_plain);
    ct.push_back(std::move(result));
    return ct;
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWEGadgetCiphertext SchemeImpl<Poly>::RLWEGadgetEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain) {
    RLWEGadgetCiphertext ct(params.g.size());
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(params.g.size()))
    for (usint i = 0; i < params.g.size(); ++i) {
        const auto &t = params.g[i];
        ct[i] = RLWEEncrypt(t * m, sk, q_plain);
    }
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
    for (auto &si : sk) {
        result.push_back(RLWEGadgetEncrypt(si, skN, params.Q));
    }
    return result;
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K) {
    usint k = ct.size() - 1;
    usint kN = K[0][0].size() - 1;
    usint l = params.g.size();
    RLWECiphertext result(kN + 1);
    for (usint i = 0; i <= kN; ++i) {
        result[i] = Poly(params.poly, EVALUATION, true);
    }
    result[kN] += ct[k];
    for (usint i = 0; i < k; ++i) {
        auto t = result[kN];
        t.SetFormat(COEFFICIENT);
        auto a = ct[i];
        std::vector<Poly> lista(l);
        a.SetFormat(COEFFICIENT);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(l))
        for (usint j = 0; j < l; ++j) {
            const auto &t = params.g[j];
            Poly a0(params.poly, COEFFICIENT, true);
            for (usint m = 0; m < a.GetLength(); ++m) {
                a0[m] = (a[m] / t) % params.Bks;
            }
            a0.SetFormat(EVALUATION);
            lista[j] = a0;
        }

        for (usint j = 0; j < l; ++j) {
            for (usint m = 0; m <= kN; ++m) {
                result[m] -= lista[j] * K[i][j][m];
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
    return std::make_pair(RLWEGadgetEncrypt(m * s, sk, params.Q), RLWEGadgetEncrypt(m, sk, params.Q));
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::Mult(Poly a, RLWEGadgetCiphertext ct) {
    usint kN = ct[0].size() - 1;
    usint l = params.g.size();
    a.SetFormat(COEFFICIENT);
    RLWECiphertext result(kN + 1);
    for (usint i = 0; i <= kN; ++i) {
        result[i] = Poly(params.poly, EVALUATION, true);
        std::vector<Poly> delta(l);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(l))
        for (usint j = 0; j < l; ++j) {
            const auto &t = params.g[j];
            Poly a0(params.poly, COEFFICIENT, true);
            for (usint m = 0; m < a.GetLength(); ++m) {
                a0[m] = (a[m] / t) % params.Bks;
            }
            a0.SetFormat(EVALUATION);
            delta[j] = a0 * ct[j][i];
        }
        for (usint j = 0; j < l; ++j) {
            result[i] += delta[j];
        }
    }
    return result;
}

template <typename Poly>
typename SchemeImpl<Poly>::RLWECiphertext SchemeImpl<Poly>::ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW) {
    if (ct.size() != 2) {
        throw std::runtime_error("ExtMult requires a ciphertext with two components");
    }
    Poly ra(params.poly, EVALUATION, true);
    Poly rb(params.poly, EVALUATION, true);
    usint l = ctGSW.first.size();
    Poly a = ct[0].Negate(), b = ct[1];
    a.SetFormat(COEFFICIENT);
    b.SetFormat(COEFFICIENT);
    std::vector<Poly> deltaa(l), deltab(l);
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(l))
    for (usint i = 0; i < l; ++i) {
        const auto &t = params.g[i];
        Poly ai(params.poly, COEFFICIENT, true);
        Poly bi(params.poly, COEFFICIENT, true);
        for (usint j = 0; j < a.GetLength(); ++j) {
            ai[j] = (a[j] / t) % params.Bks;
            bi[j] = (b[j] / t) % params.Bks;
        }
        ai.SetFormat(EVALUATION);
        bi.SetFormat(EVALUATION);
        deltaa[i] = ai * ctGSW.first[i][0] + bi * ctGSW.second[i][0];
        deltab[i] = ai * ctGSW.first[i][1] + bi * ctGSW.second[i][1];
    }
    for (usint i = 0; i < l; i++) {
        ra += deltaa[i];
        rb += deltab[i];
    }
    return {std::move(ra), std::move(rb)};
}

#endif