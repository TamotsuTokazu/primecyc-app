#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

#include <chrono>

#define START_TIMER start = std::chrono::system_clock::now()
#define END_TIMER std::cout << "Time: " << (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()) << std::endl

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

using lbcrypto::BigInteger;

namespace par {

const Integer t("64");
const usint n = 100;
const usint p0 = 1153;
const usint p1 = 1297;
const usint pq = p0 * p1;
const usint Bks = 1 << 4;
const Integer Q("1152920977604149249");

const Integer rootOfUnity("739447795444923848");
const Integer rootOfUnity0("587824393635068054");
const Integer rootOfUnity1("950015890856494149");

const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, 0, 0);
const auto pp1 = std::make_shared<ILParams>(p1, Q, rootOfUnity1, 0, 0);
const auto ppq = std::make_shared<ILParams>(pq, Q, rootOfUnity, 0, 0);

}

Poly Tensor(const Poly &a, const Poly &b);
RLWEKey TensorKey(const RLWEKey &skp, const RLWEKey &skq);
RLWECiphertext TensorCt(const RLWECiphertext &ap, const RLWECiphertext &aq);

Vector TracePqToP(const Vector &a, usint p, usint q, Integer Q);
Poly TracePqToP(const Poly &a);
Integer TracePtoZ(const Poly &a);

Poly ConstructF(const std::vector<usint> &f);
Poly ConstructFP(const std::vector<usint> &f, std::shared_ptr<ILParams> params);
Poly ConstructKey(Vector sk, std::shared_ptr<ILParams> params);
Vector RetrieveLWECt(const RLWECiphertext &ct);

BigInteger Resize(BigInteger x, BigInteger q, BigInteger Q);

int main() {

    std::chrono::time_point<std::chrono::system_clock> start;

    primecyc::RaderFFTNat<Vector>::m_enabled[par::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p1] = true;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(par::pq);
    Vector sk = dugpq.GenerateVector(par::n);
    Vector a = dugpq.GenerateVector(par::n);
    usint m_plain = 1;
    Integer b = m_plain;
    std::vector<usint> f_plain(par::pq, 0);
    f_plain[1] = 1;

    std::vector<usint> f_arr(par::pq, 0);
    for (usint i = 0; i < par::pq; i++) {
        f_arr[i] = f_plain[(par::pq - i) % par::pq];
    }

    for (usint i = 0; i < par::n; i++) {
        b += a[i] * sk[i];
    }

    std::cout << "Stage 1: Prime KeyGen" << std::endl;

    START_TIMER;

    Scheme sc0{{par::pp0, par::p0, par::Q, par::Bks}, sk};
    Scheme sc1{{par::pp1, par::p1, par::Q, par::Bks}, sk};

    END_TIMER;

    std::cout << "Stage 2: Prime Eval" << std::endl;

    START_TIMER;

    auto mult = Integer(par::pq).ModInverse(par::t);

    auto ct0 = sc0.Process(a, b, par::t, mult);
    auto ct1 = sc1.Process(a, b, par::t, mult);

    primecyc::TensorFFTNat<Vector>::m_factor[par::pq] = {par::p0, par::p1};

    END_TIMER;

    std::cout << "Stage 3: Tensor" << std::endl;

    START_TIMER;

    auto ct = TensorCt(ct0, ct1);
    auto skk = TensorKey({sc0.skp}, {sc1.skp});

    Scheme scheme_tensor({par::ppq, par::pq, par::Q, par::Bks});

    Poly skpq = Poly(par::ppq, COEFFICIENT, true);
    for (usint i = 0; i < par::n; i++) {
        if (i == 0) {
            for (usint j = 0; j < skpq.GetLength(); j++) {
                skpq[j] = sk[i];
            }
        } else {
            usint t = (i * par::p1 % par::p0 - 1) * (par::p1 - 1);
            for (usint k = 1; k < par::p1; k++) {
                skpq[t + k - 1].ModSubEq(sk[i], par::Q);
            }
        }
    }

    END_TIMER;

    std::cout << "Stage 4: Tensor Switching KeyGen" << std::endl;

    START_TIMER;

    skpq.SetFormat(EVALUATION);
    auto tensor_ksk = scheme_tensor.KeySwitchGen(skk, {skpq});

    END_TIMER;

    std::cout << "Stage 5: Tensored Key Switching" << std::endl;

    START_TIMER;

    auto tensor_ct = scheme_tensor.KeySwitch(ct, tensor_ksk);

    auto f = ConstructF(f_arr);
    f.SetFormat(EVALUATION);
    tensor_ct[0] *= f;
    tensor_ct[1] *= f;

    END_TIMER;

    std::cout << "Stage 6: Trace" << std::endl;

    START_TIMER;

    auto cta = tensor_ct[0];
    auto ctb = tensor_ct[1];
    cta.SetFormat(COEFFICIENT);
    ctb.SetFormat(COEFFICIENT);
    cta = TracePqToP(cta);
    ctb = TracePqToP(ctb);

    auto skp = ConstructKey(sk, par::pp0);

    cta.SetFormat(EVALUATION);
    ctb.SetFormat(EVALUATION);
    skp.SetFormat(EVALUATION);

    auto fp = ConstructFP(f_arr, par::pp0);
    fp.SetFormat(EVALUATION);
    auto inv0 = Integer(par::p1).ModInverse(par::p0).ConvertToInt();
    auto gal_t0 = sc0.GaloisConjugate(sc0.skp, inv0);
    ct0[0] *= fp;
    ct0[1] *= fp;
    auto ksk_extra0 = sc0.KeySwitchGen({gal_t0}, {skp});
    auto ct_extra0 = sc0.GaloisConjugate(ct0, inv0);
    ct_extra0 = sc0.KeySwitch(ct_extra0, ksk_extra0);

    cta += ct_extra0[0];
    ctb += ct_extra0[1];

    auto skq = ConstructKey(sk, par::pp1);
    skq.SetFormat(EVALUATION);
    auto fq = ConstructFP(f_arr, par::pp1);
    fq.SetFormat(EVALUATION);
    auto inv1 = Integer(par::p0).ModInverse(par::p1).ConvertToInt();
    auto gal_t1 = sc1.GaloisConjugate(sc1.skp, inv1);
    ct1[0] *= fq;
    ct1[1] *= fq;
    auto ksk_extra1 = sc1.KeySwitchGen({gal_t1}, {skq});
    auto ct_extra1 = sc1.GaloisConjugate(ct1, inv1);
    ct_extra1 = sc1.KeySwitch(ct_extra1, ksk_extra1);

    auto plain = sc0.RLWEDecrypt({cta, ctb}, {skp}, par::Q);
    sc0.ModSwitch(plain, par::t);

    cta.SetFormat(COEFFICIENT);
    ctb.SetFormat(COEFFICIENT);

    ct_extra1[0].SetFormat(COEFFICIENT);
    ct_extra1[1].SetFormat(COEFFICIENT);

    auto lwect = RetrieveLWECt({cta, ctb}) + RetrieveLWECt(ct_extra1);
    lwect[par::n].ModAddEq((par::Q / par::t).ModMul(mult, par::Q), par::Q);

    END_TIMER;

    Integer result = lwect[par::n];
    for (usint i = 0; i < par::n; i++) {
        result.ModSubEq(lwect[i].ModMulEq(sk[i], par::Q), par::Q);
    }

    std::cout << "decrypted: " << Resize(result, par::t, par::Q) << std::endl;
    std::cout << "expected: " << f_plain[m_plain] << std::endl;


    return 0;
}

Vector TracePqToP(const Vector &a, usint p, usint q, Integer Q) {
    Vector b(p - 1, Q);
    for (usint i = 1; i < p; i++) {
        usint t = ((q * i) % p - 1) * (q - 1);
        for (usint j = 1; j < q; j++) {
            b[i - 1].ModSubEq(a[t + j - 1], Q);
        }
    }
    return b;
}

Poly TracePqToP(const Poly &a) {
    Vector b = TracePqToP(a.GetValues(), par::p0, par::p1, par::Q);
    Poly c(par::pp0, COEFFICIENT, true);
    for (usint i = 1; i < par::p0 - 1; i++) {
        c[i] = b[i - 1].ModSub(b[par::p0 - 2], par::Q);
    }
    c[0] = par::Q - b[par::p0 - 2];
    return c;
}

Poly Tensor(const Poly &a, const Poly &b) {
    Poly c(par::ppq, EVALUATION, true);
    for (usint i = 1; i < par::p0; i++) {
        for (usint j = 1; j < par::p1; j++) {
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
    for (usint i = 0; i < par::p0 - 1; i++) {
        p1[i] = 1;
    }
    for (usint i = 0; i < par::p1 - 1; i++) {
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

Integer TracePtoZ(const Poly &a) {
    auto p = a.GetCyclotomicOrder();
    Integer c = a[0].ModMul(p - 1, par::Q);
    for (usint i = 1; i < p - 1; i++) {
        c.ModSubEq(a[i], par::Q);
    }
    return c;
}

Poly ConstructF(const std::vector<usint> &f) {
    Poly result(par::ppq, COEFFICIENT, true);
    for (usint k = 0; k < result.GetLength(); k++) {
        result[k] = f[0];
    }
    for (usint k = 1; k < f.size(); k++) {
        usint i = k % par::p0;
        usint j = k % par::p1;
        if (i == 0) {
            for (usint l = 1; l < par::p0; l++) {
                result[(l - 1) * (par::p1 - 1) + j - 1].ModSubEq(f[k], par::Q);
            }
        } else if (j == 0) {
            for (usint l = 1; l < par::p1; l++) {
                result[(i - 1) * (par::p1 - 1) + l - 1].ModSubEq(f[k], par::Q);
            }
        } else {
            result[(i - 1) * (par::p1 - 1) + j - 1].ModAddEq(f[k], par::Q);
        }
    }
    return result;
}

Poly ConstructFP(const std::vector<usint> &f, std::shared_ptr<ILParams> params) {
    Poly result(params, COEFFICIENT, true);
    auto p = params->GetCyclotomicOrder();
    Integer t;
    for (usint k = 0; k < f.size(); k++) {
        if (k % p == p - 1) {
            t.ModAddEq(f[k], par::Q);
        } else {
            result[k % p].ModAddEq(f[k], par::Q);
        }
    }
    for (usint k = 0; k < p - 1; k++) {
        result[k].ModSubEq(t, par::Q);
    }
    return result;
}

Poly ConstructKey(Vector sk, std::shared_ptr<ILParams> params) {
    Poly result(params, COEFFICIENT, true);
    for (usint i = 0; i < par::n; i++) {
        result[i] = sk[i];
    }
    return result;
}

Vector RetrieveLWECt(const RLWECiphertext &ct) {
    Vector result(par::n + 1);
    Integer a_sum = 0;
    auto &cta = ct[0];
    usint p = cta.GetCyclotomicOrder();

    for (usint i = 0; i < cta.GetLength(); i++) {
        a_sum.ModAddEq(cta[i], par::Q);
    }
    for (usint i = 0; i < par::n; i++) {
        if (i == 0) {
            result[i] = cta[0].ModMul(p, par::Q).ModSub(a_sum, par::Q);
        } else if (i == 1) {
            result[i] = Integer(0).ModSub(a_sum, par::Q);
        } else {
            result[i] = cta[p - i].ModMul(p, par::Q).ModSub(a_sum, par::Q);
        }
    }

    result[par::n] = TracePtoZ(ct[1]);
    return result;
}

BigInteger Resize(BigInteger x, BigInteger q, BigInteger Q) {
    BigInteger halfQ = Q / 2;
    x = (x * q + halfQ) / Q;
    if (x >= q) {
        x -= q;
    }
    return x;
}