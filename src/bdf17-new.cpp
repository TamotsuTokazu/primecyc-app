#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

using Poly = lbcrypto::NativePoly;
using Scheme = SchemeImpl<Poly>;
using Vector = Scheme::Vector;
using Integer = Scheme::Integer;
using ILParams = Scheme::ILParams;
using Params = Scheme::Params;

using RLWECiphertext = Scheme::RLWECiphertext;
using RLWEKey = Scheme::RLWEKey;
using RLWESwitchingKey = Scheme::RLWESwitchingKey;
using RGSWCiphertext = Scheme::RGSWCiphertext;

namespace par {
    const Integer t = 1 << 6;
    const uint32_t n = 100;
    const uint32_t p0 = 1153;
    const uint32_t p1 = 1297;
    const uint32_t pq = p0 * p1;
    const uint32_t Bks = 1 << 6;
    
    const Integer Q = 36066736134770689;
    const Integer rootOfUnity = 4364918564594134;
    const Integer rootOfUnity0 = 28679241126083710;
    const Integer rootOfUnity1 = 9730598941305417;

    const Integer bigModulus = Q;
    const Integer bigRootOfUnity0 = rootOfUnity0;
    const Integer bigRootOfUnity1 = rootOfUnity1;

    const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, bigModulus, bigRootOfUnity0);
    const auto pp1 = std::make_shared<ILParams>(p1, Q, rootOfUnity1, bigModulus, bigRootOfUnity1);
    const auto ppq = std::make_shared<ILParams>(pq, Q, rootOfUnity, bigModulus, rootOfUnity);

}

Poly Tensor(const Poly &a, const Poly &b);
RLWEKey TensorKey(const RLWEKey &skp, const RLWEKey &skq);
RLWECiphertext TensorCt(const RLWECiphertext &ap, const RLWECiphertext &aq);

Vector TracePqToP(const Vector &a, uint32_t p, uint32_t q, Integer Q);
Poly TracePqToP(const Poly &a);

int main() {
    uint32_t pq = par::p0 * par::p1;

    primecyc::RaderFFTNat<Vector>::m_enabled[par::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p1] = true;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(pq);
    Vector sk = dugpq.GenerateVector(par::n);
    Vector a = dugpq.GenerateVector(par::n);
    Integer b = 0;

    for (uint32_t i = 0; i < par::n; i++) {
        b += a[i] * sk[i];
    }

    std::cout << "Stage 1" << std::endl;

    Scheme sc0{{par::pp0, par::p0, par::Q, par::Bks}, sk};
    Scheme sc1{{par::pp1, par::p1, par::Q, par::Bks}, sk};

    std::cout << "Stage 2" << std::endl;

    auto ct0 = sc0.Process(a, b, par::t);
    auto ct1 = sc1.Process(a, b, par::t);

    primecyc::TensorFFTNat<Vector>::m_factor[pq] = {par::p0, par::p1};

    std::cout << "Stage 3" << std::endl;

    auto ct = TensorCt(ct0, ct1);
    auto skk = TensorKey({sc0.skp}, {sc1.skp});

    std::cout << "Stage 4" << std::endl;

    Scheme scheme_tensor({par::ppq, pq, par::Q, par::Bks});

    Poly skpq = Poly(par::ppq, COEFFICIENT, true);
    for (uint32_t i = 0; i < par::n; i++) {
        if (i == 0) {
            for (uint32_t j = 0; j < skpq.GetLength(); j++) {
                skpq[j] = sk[i];
            }
        }
        u_int32_t t = i * par::p1 % par::p0;
        for (u_int32_t k = 1; k < par::p1; k++) {
            skpq[t * (par::p1 - 1) + k - 1].ModSub(sk[i], par::Q);
        }
    }
    skpq.SetFormat(EVALUATION);

    std::cout << "Stage 5" << std::endl;

    auto tensor_ksk = scheme_tensor.KeySwitchGen(skk, {skpq});

    std::cout << "Stage 6" << std::endl;

    auto tensor_ct = scheme_tensor.KeySwitch(ct, tensor_ksk);

    std::cout << "Stage 7" << std::endl;

    skpq.SetFormat(COEFFICIENT);
    auto cta = tensor_ct[0];
    auto ctb = tensor_ct[1];
    cta.SetFormat(COEFFICIENT);
    ctb.SetFormat(COEFFICIENT);
    cta = TracePqToP(cta);
    ctb = TracePqToP(ctb);
    skpq = TracePqToP(skpq) * Integer(par::p1 - 1).ModInverse(par::Q);
    cta.SetFormat(EVALUATION);
    ctb.SetFormat(EVALUATION);
    skpq.SetFormat(EVALUATION);

    auto plain = sc0.RLWEDecrypt({cta, ctb}, {skpq}, par::Q);
    sc0.ModSwitch(plain, par::t);
    std::cout << "decrypted: " << plain << std::endl;

    return 0;
}

Vector TracePqToP(const Vector &a, uint32_t p, uint32_t q, Integer Q) {
    Vector b(p - 1, Q);
    for (uint32_t i = 1; i < p; i++) {
        uint32_t t = ((q * i) % p - 1) * (q - 1);
        for (uint32_t j = 1; j < q; j++) {
            b[i - 1].ModSubEq(a[t + j - 1], Q);
        }
    }
    return b;
}

Poly TracePqToP(const Poly &a) {
    Vector b = TracePqToP(a.GetValues(), par::p0, par::p1, par::Q);
    Poly c(par::pp0, COEFFICIENT, true);
    for (uint32_t i = 1; i < par::p0 - 1; i++) {
        c[i] = b[i - 1].ModSub(b[par::p0 - 2], par::Q);
    }
    c[0] = par::Q - b[par::p0 - 2];
    return c;
}

Poly Tensor(const Poly &a, const Poly &b) {
    Poly c(par::ppq, EVALUATION, true);
    for (uint32_t i = 1; i < par::p0; i++) {
        for (uint32_t j = 1; j < par::p1; j++) {
            c[(i - 1) * (par::p1 - 1) + j - 1] = a[i - 1].ModMul(b[j - 1], par::Q);
        }
    }
    return c;
}

RLWEKey TensorKey(const RLWEKey &skp, const RLWEKey &skq) {
    auto skp0 = skp[0];
    auto skq0 = skq[0];
    Poly p1(par::pp0, EVALUATION, true);
    Poly q1(par::pp1, EVALUATION, true);
    for (uint32_t i = 0; i < par::p0 - 1; i++) {
        p1[i] = 1;
    }
    for (uint32_t i = 0; i < par::p1 - 1; i++) {
        q1[i] = 1;
    }
    RLWEKey sk;
    sk.push_back(Tensor(skp0, skq0).Negate());
    sk.push_back(Tensor(skp0, q1));
    sk.push_back(Tensor(p1, skq0));
    return sk;
}

RLWECiphertext TensorCt(const RLWECiphertext &ap, const RLWECiphertext &aq) {
    RLWECiphertext c;
    Integer z = par::Q - par::t;
    c.push_back(z * Tensor(ap[0], aq[0]));
    c.push_back(z * Tensor(ap[0], aq[1]));
    c.push_back(z * Tensor(ap[1], aq[0]));
    c.push_back(z * Tensor(ap[1], aq[1]));
    return c;
}