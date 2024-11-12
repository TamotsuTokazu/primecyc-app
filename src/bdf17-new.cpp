#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

#include <chrono>

#define START_TIMER start = std::chrono::system_clock::now()
#define END_TIMER std::cout << "Time: " << (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()) << std::endl

using Poly = lbcrypto::NativePoly;
using Scheme = BDF17SchemeImpl<Poly>;
using Vector = Scheme::Vector;
using Integer = Scheme::Integer;
using ILParams = Scheme::ILParams;

using RLWECiphertext = Scheme::RLWECiphertext;
using RLWEKey = Scheme::RLWEKey;
using RLWESwitchingKey = Scheme::RLWESwitchingKey;
using RGSWCiphertext = Scheme::RGSWCiphertext;

using lbcrypto::BigInteger;

namespace par {

const Integer t("8 ");
const usint n = 600;
const usint p0 = 1153;
const usint p1 = 1297;
const usint pq = p0 * p1;
const usint Bks = 1 << 8;
const Integer Q("576485048298018817");

const Integer rootOfUnity("537122826862086034");
const Integer rootOfUnity0("397835389725933388");
const Integer rootOfUnity1("440043149905686213");

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
Vector RetrieveLWECt(const RLWECiphertext &ct, usint len);

BigInteger Resize(BigInteger x, BigInteger q, BigInteger Q);

int main() {

    std::chrono::time_point<std::chrono::system_clock> start;

    primecyc::RaderFFTNat<Vector>::m_enabled[par::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p1] = true;
    primecyc::TensorFFTNat<Vector>::m_factor[par::pq] = {par::p0, par::p1};

    Poly x10 = Poly(par::pp0, COEFFICIENT, true);
    Poly x11 = Poly(par::pp1, COEFFICIENT, true);
    x10[0] = par::Q - 1;
    x11[0] = par::Q - 1;
    x10[1] = 1;
    x11[1] = 1;
    x10.SetFormat(EVALUATION);
    x11.SetFormat(EVALUATION);

    Poly one0 = Poly(par::pp0, COEFFICIENT, true);
    Poly one1 = Poly(par::pp1, COEFFICIENT, true);
    one0[0] = 1;
    one1[0] = 1;
    one0.SetFormat(EVALUATION);
    one1.SetFormat(EVALUATION);

    Poly x1pq = Tensor(x10, x11);

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(par::pq);

    std::cout << "Stage 0: KeyGen" << std::endl;

    START_TIMER;

    Vector sk = dugpq.GenerateVector(par::n);
    Scheme sc0{{par::pp0, par::p0, par::Q, par::Bks}, x10, sk};
    Scheme sc1{{par::pp1, par::p1, par::Q, par::Bks}, x11, sk};

    Scheme scheme_tensor({par::ppq, par::pq, par::Q, par::Bks}, x1pq);
    auto skk = TensorKey({sc0.skp}, {sc1.skp});
    auto sk_tensor_0 = Tensor(sc0.skp, one1);
    auto tensor_ksk0 = scheme_tensor.KeySwitchGen(skk, {sk_tensor_0});

    END_TIMER;

    for (usint numRound = 0; numRound < 10; numRound++) {

        Vector a = dugpq.GenerateVector(par::n);
        usint m_plain = numRound;
        Integer b = m_plain * (usint)(0.5 + par::pq / par::t.ConvertToDouble());
        std::vector<usint> f_plain(par::t.ConvertToInt(), 0);
        for (usint i = 0; i < par::t.ConvertToInt(); i++) {
            f_plain[i] = i & 1;
        }

        std::vector<usint> f_ct(par::pq, 0);
        for (usint i = 0; i < par::pq; i++) {
            f_ct[i] = f_plain[(usint)(0.5 + par::t.ConvertToDouble() * i / par::pq) % par::t.ConvertToInt()];
        }
        for (usint i = 1, j = par::pq - 1; i < j; i++, j--) {
            std::swap(f_ct[i], f_ct[j]);
        }

        for (usint i = 0; i < par::n; i++) {
            b += a[i] * sk[i];
        }

        std::cout << "Stage 1: Prime Eval" << std::endl;

        START_TIMER;

        auto ct0 = sc0.Process(a, b, par::t);
        auto ct1 = sc1.Process(a, b, par::t);

        END_TIMER;

        std::cout << "Stage 2: Tensor" << std::endl;

        START_TIMER;

        auto ct = TensorCt(ct0, ct1);

        END_TIMER;

        std::cout << "Stage 3: Tensored Key Switching" << std::endl;

        START_TIMER;

        auto tensor_ct = scheme_tensor.KeySwitch(ct, tensor_ksk0);

        END_TIMER;

        std::cout << "Stage 4: Function Extraction" << std::endl;

        START_TIMER;

        auto f = ConstructF(f_ct);
        f.SetFormat(EVALUATION);
        tensor_ct[0] *= f;
        tensor_ct[1] *= f;

        auto cta = tensor_ct[0];
        auto ctb = tensor_ct[1];
        cta.SetFormat(COEFFICIENT);
        ctb.SetFormat(COEFFICIENT);

        cta = TracePqToP(cta);
        ctb = TracePqToP(ctb);

        cta.SetFormat(EVALUATION);
        ctb.SetFormat(EVALUATION);

        cta = sc0.GaloisConjugate(cta, par::p1 % par::p0);
        ctb = sc0.GaloisConjugate(ctb, par::p1 % par::p0);

        auto fp = ConstructFP(f_ct, par::pp0);
        fp.SetFormat(EVALUATION);

        cta += ct0[0] * fp;
        ctb += ct0[1] * fp;

        sc0.skp.SetFormat(COEFFICIENT);
        cta.SetFormat(COEFFICIENT);
        ctb.SetFormat(COEFFICIENT);

        auto fq = ConstructFP(f_ct, par::pp1);
        fq.SetFormat(EVALUATION);
        ct1[0] *= fq;
        ct1[1] *= fq;
        ct1[0].SetFormat(COEFFICIENT);
        ct1[1].SetFormat(COEFFICIENT);

        auto lwect = RetrieveLWECt({cta, ctb}, par::n) + RetrieveLWECt(ct1, par::n);

        lwect[par::n] += par::Q / par::t;

        END_TIMER;

        Integer result = lwect[par::n];
        for (usint i = 0; i < par::n; i++) {
            result.ModSubEq(lwect[i].ModMulEq(sk[i], par::Q), par::Q);
        }
        result.ModMulEq(Integer(par::pq).ModInverse(par::Q), par::Q);

        std::cout << "decrypted: " << Resize(result, par::t, par::Q) << std::endl;
        std::cout << "expected: " << f_plain[m_plain] << std::endl;

    }

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

Vector RetrieveLWECt(const RLWECiphertext &ct, usint len) {
    Vector result(len + 1);
    Integer a_sum = 0;
    auto &cta = ct[0];
    usint p = cta.GetCyclotomicOrder();

    for (usint i = 0; i < cta.GetLength(); i++) {
        a_sum.ModAddEq(cta[i], par::Q);
    }

    for (usint i = 0; i < len; i++) {
        if (i == 0) {
            result[i] = cta[0].ModMul(p, par::Q).ModSub(a_sum, par::Q);
        } else if (i == 1) {
            result[i] = Integer(0).ModSub(a_sum, par::Q);
        } else {
            result[i] = cta[p - i].ModMul(p, par::Q).ModSub(a_sum, par::Q);
        }
    }

    result[len] = TracePtoZ(ct[1]);
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