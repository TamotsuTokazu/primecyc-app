#ifndef _BDF17SCHEME_IMPL_H_
#define _BDF17SCHEME_IMPL_H_

template <typename Poly>
BDF17SchemeImpl<Poly>::BDF17SchemeImpl(Params p, Poly x1_, Vector sk) : SchemeImpl<Poly>(p, false), bk(sk.GetLength()), x1(p.poly, EVALUATION, false) {
    x1 = x1_;
    sk.SetModulus(p.p);
    sk.ModEq(p.p);
    if (sk.GetLength() != 0) {
        this->skp = Poly(p.poly, COEFFICIENT, true);
        for (usint i = 0; i < sk.GetLength(); i++) {
            this->skp[i] = sk[i];
        }
        this->GaloisKeyGen();
#pragma omp parallel for num_threads(lbcrypto::OpenFHEParallelControls.GetThreadLimit(sk.GetLength()))
        for (usint i = 0; i < sk.GetLength(); i++) {
            Poly m(p.poly, COEFFICIENT, true);
            auto t = sk[i].ConvertToInt();
            if (t < p.p - 1) {
                m[sk[i].ConvertToInt()] = 1;
            } else {
                for (usint j = 0; j < p.p - 1; j++) {
                    m[j] = Integer(p.poly->GetModulus() - 1);
                }
            }
            m.SetFormat(EVALUATION);
            bk[i] = this->RGSWEncrypt(m, {this->skp});
        }
    }
}

template <typename Poly>
typename BDF17SchemeImpl<Poly>::RLWECiphertext BDF17SchemeImpl<Poly>::RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain) {
    const auto &params = this->params;

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
    result += e * x1 + m * (params.Q / q_plain);
    ct.push_back(std::move(result));
    return ct;
}

template <typename Poly>
typename BDF17SchemeImpl<Poly>::RLWECiphertext BDF17SchemeImpl<Poly>::Process(Vector a, Integer b, Integer q_plain) {
    const auto &params = this->params;
    a.SetModulus(params.p);
    a.ModEq(params.p);
    b.ModEq(params.p);
    Poly ca(params.poly, EVALUATION, true);
    Poly cb(params.poly, COEFFICIENT, true);
    if (b == params.p - 1) {
        for (usint i = 0; i < params.p - 1; i++) {
            cb[i] = params.Q - params.Q / q_plain;
        }
    } else {
        cb[b.ConvertToInt()] = params.Q / q_plain;
    }
    cb.SetFormat(EVALUATION);
    RLWECiphertext c{ca, cb};
    Integer t = 1;
    for (usint i = 0; i < bk.size(); i++) {
        a[i] = params.p - a[i];
        if (a[i] != params.p) {
            t.ModMulEq(a[i].ModInverse(params.p), params.p);
            if (t != 1) {
                c = this->GaloisConjugate(c, t.ConvertToInt());
                c = this->KeySwitch(c, this->ksk_galois[t.ConvertToInt()]);
            }
            c = this->ExtMult(c, bk[i]);
            t = a[i];
        }
    }
    if (t != 1) {
        c = this->GaloisConjugate(c, t.ConvertToInt());
        c = this->KeySwitch(c, this->ksk_galois[t.ConvertToInt()]);
    }
    return c;
}

#endif // _BDF17SCHEME_IMPL_H_