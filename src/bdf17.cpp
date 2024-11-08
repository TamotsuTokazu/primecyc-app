#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

#include <chrono>

#define START_TIMER start = std::chrono::system_clock::now()
#define END_TIMER std::cout << "Time: " << (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()) << std::endl

using Poly = lbcrypto::Poly;
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
    const Integer t("65");
    const usint n = 600;
    const usint p0 = 1153;
    const usint p1 = 1297;
    const usint pq = p0 * p1;
    const usint Bks = 1 << 6;
    const Integer Q("72057595543256441");


    const Integer rootOfUnitypq("67697845297299456");
    const Integer rootOfUnity0("4261944152330832");
    const Integer rootOfUnity1("45688242756488110");


    const Integer bigModulus("22300745198530623141535718272648361904439297");
    const Integer bigRootOfUnity0("10582516286048940371395447171034973190070517");
    const Integer bigRootOfUnity1("10582516286048940371395447171034973190070517");
    const Integer bigRootOfUnitypq("1773029494350011920780512139097380038959158");
    // const Integer Q = lbcrypto::LastPrime<Integer>(56, pq * t.ConvertToInt());
    // const Integer rootOfUnitypq = lbcrypto::RootOfUnity<Integer>(2 * pq, Q);
    // const Integer rootOfUnity0 = rootOfUnitypq.ModExp(p1, Q);
    // const Integer rootOfUnity1 = rootOfUnitypq.ModExp(p0, Q);

    // const usint nttSize0 = 1 << (int)std::ceil(std::log2(2 * p0 - 1));
    // const usint nttSize1 = 1 << (int)std::ceil(std::log2(2 * p1 - 1));
    // const usint nttSizepq = 1 << (int)std::ceil(std::log2(2 * pq - 1));

    // const Integer bigModulus = lbcrypto::FirstPrime<Integer>(1 + (int)(std::ceil(std::log2(2 * pq - 1) + 2 * std::log2(Q.ConvertToInt()))), nttSizepq);
    // const Integer bigRootOfUnity0 = lbcrypto::RootOfUnity<Integer>(nttSize0, bigModulus);
    // const Integer bigRootOfUnity1 = lbcrypto::RootOfUnity<Integer>(nttSize1, bigModulus);
    // const Integer bigRootOfUnitypq = lbcrypto::RootOfUnity<Integer>(nttSizepq, bigModulus);

    const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, bigModulus, bigRootOfUnity0);
    const auto pp1 = std::make_shared<ILParams>(p1, Q, rootOfUnity1, bigModulus, bigRootOfUnity1);
    const auto ppq = std::make_shared<ILParams>(pq, Q, rootOfUnitypq, bigModulus, bigRootOfUnitypq);

}


Poly Tensor(const Poly &a, const Poly &b);
RLWEKey TensorKey(const RLWEKey &skp, const RLWEKey &skq);
RLWECiphertext TensorCt(const RLWECiphertext &ap, const RLWECiphertext &aq);

Poly TracePqToP(const Poly &a);
Integer TracePtoZ(const Poly &a);

int main() {

    std::chrono::time_point<std::chrono::system_clock> start;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(par::pq);
    Vector sk = dugpq.GenerateVector(par::n);
    Vector a = dugpq.GenerateVector(par::n);
    Integer b = 0;
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

    auto ct0 = sc0.Process(a, b, par::t);
    auto ct1 = sc1.Process(a, b, par::t);

    END_TIMER;

    std::cout << "Stage 3: Tensor" << std::endl;

    START_TIMER;

    auto ct = TensorCt(ct0, ct1);
    auto skk = TensorKey({sc0.skp}, {sc1.skp});

    Scheme scheme_tensor({par::ppq, par::pq, par::Q, par::Bks});
    Poly skpq = Poly(par::ppq, COEFFICIENT, true);
    for (usint i = 0; i < par::n; i++) {
        skpq[i * par::p1] = sk[i];
    }
    skpq.SetFormat(EVALUATION);

    END_TIMER;

    std::cout << "Stage 4: Tensor Switching KeyGen" << std::endl;

    START_TIMER;

    auto tensor_ksk = scheme_tensor.KeySwitchGen(skk, {skpq});

    END_TIMER;

    std::cout << "Stage 5: Tensored Key Switching" << std::endl;

    START_TIMER;

    auto tensor_ct = scheme_tensor.KeySwitch(ct, tensor_ksk);

    END_TIMER;

    std::cout << "Stage 6: Trace" << std::endl;

    START_TIMER;

    auto cta = tensor_ct[0];
    auto ctb = tensor_ct[1];
    cta.SetFormat(COEFFICIENT);
    ctb.SetFormat(COEFFICIENT);
    skpq.SetFormat(COEFFICIENT);
    cta = TracePqToP(cta);
    ctb = TracePqToP(ctb);
    skpq = TracePqToP(skpq) * Integer(par::p1 - 1).ModInverse(par::Q);
    cta.SetFormat(EVALUATION);
    ctb.SetFormat(EVALUATION);
    skpq.SetFormat(EVALUATION);

    auto plain = scheme_tensor.RLWEDecrypt({cta, ctb}, {skpq}, par::Q);
    scheme_tensor.ModSwitch(plain, par::t);

    END_TIMER;

    std::cout << "decrypted: " << plain << std::endl;

    Integer z = TracePtoZ(plain);
    z = (z * par::t + (par::Q / 2)) / par::Q;

    std::cout << "decrypted: " << z << std::endl;

    return 0;
}

Poly Tensor(const Poly &a, const Poly &b) {
    // assume a follows p0 and b follows p1
    Poly c(par::ppq, EVALUATION, true);
    std::vector<usint> idx(par::pq);
    usint cur = 0;
    for (usint i = 1; i < par::pq; i++) {
        if (i % par::p0 != 0 && i % par::p1 != 0) {
            idx[i] = cur++;
        }
    }
    for (usint i = 1; i < par::p0; i++) {
        for (usint j = 1; j < par::p1; j++) {
            c[idx[(i * par::p1 + j * par::p0) % par::pq]] = a[i - 1].ModMul(b[j - 1], par::Q);
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

Poly TracePqToP(const Poly &a) {
    usint n = a.GetLength();
    Poly c(par::pp0, COEFFICIENT, true);
    for (usint i = 0; i < par::p0; i++) {
        for (usint j = 0; j < par::p1; j++) {
            auto t = (par::p1 * i + par::p0 * j) % par::pq;
            if (t < n) {
                if (j == 0) {
                    auto x = a[t].ModMul(par::p1 - 1, par::Q);
                    if (i == par::p0 - 1) {
                        for (usint k = 0; k < par::p0 - 1; k++) {
                            c[k].ModSubEq(x, par::Q);
                        }
                    } else {
                        c[i].ModAddEq(x, par::Q);
                    }
                } else {
                    if (i == par::p0 - 1) {
                        for (usint k = 0; k < par::p0 - 1; k++) {
                            c[k].ModAddEq(a[t], par::Q);
                        }
                    } else {
                        c[i].ModSubEq(a[t], par::Q);
                    }
                }
            }
        }
    }
    return c;
}

Integer TracePtoZ(const Poly &a) {
    Integer c = a[0].ModMul(par::p0 - 1, par::Q);
    for (usint i = 1; i < par::p0 - 1; i++) {
        c.ModSubEq(a[i], par::Q);
    }
    return c;
}
