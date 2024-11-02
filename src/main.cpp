#include "binfhecontext.h"
#include "params.h"
#include "rlwe.h"
#include "math/math-hal.h"

int main() {
    uint32_t pq = p::p0 * p::p1;

    primecyc::RaderFFTNat<Vector>::m_enabled[p::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[p::p1] = true;

    // print params
    std::cout << "p0: " << p::p0 << std::endl;
    std::cout << "p1: " << p::p1 << std::endl;
    std::cout << "pq: " << pq << std::endl;
    std::cout << "Q: " << p::Q << std::endl;
    std::cout << "rootOfUnity0: " << p::rootOfUnity0 << std::endl;
    std::cout << "rootOfUnity0: " << p::rootOfUnity0.ModExp(p::p0 - 1, p::Q) << std::endl;
    std::cout << "rootOfUnity0: " << p::rootOfUnity0.ModExp(p::p0, p::Q) << std::endl;
    std::cout << "rootOfUnity0: " << p::rootOfUnity0.ModExp(p::p0 * (p::p0 - 1), p::Q) << std::endl;
    std::cout << "rootOfUnity1: " << p::rootOfUnity1 << std::endl;
    std::cout << "rootOfUnity1: " << p::rootOfUnity1.ModExp(p::p1 - 1, p::Q) << std::endl;
    std::cout << "rootOfUnity1: " << p::rootOfUnity1.ModExp(p::p1, p::Q) << std::endl;
    std::cout << "rootOfUnity1: " << p::rootOfUnity1.ModExp(p::p1 * (p::p1 - 1), p::Q) << std::endl;
    std::cout << "rootOfUnitypq: " << p::rootOfUnitypq << std::endl;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(pq);
    Vector sk = dugpq.GenerateVector(p::n);
    Vector a = dugpq.GenerateVector(p::n);

    Integer b = 2;
    for (uint32_t i = 0; i < p::n; i++) {
        b += a[i] * sk[i];
    }

    Scheme sc0{{p::pp0, p::p0}, sk};
    auto ct0 = sc0.Process(a, b, p::t);
    std::cout << "ct0: " << sc0.RLWEDecrypt(ct0, {sc0.skp}, p::t) << std::endl;


    Scheme sc1{{p::pp1, p::p1}, sk};
    auto ct1 = sc1.Process(a, b, p::t);
    std::cout << "ct1: " << sc1.RLWEDecrypt(ct1, {sc1.skp}, p::t) << std::endl;

    // ChineseRemainderTransformArb<Vector>().SetCylotomicPolynomial(lbcrypto::GetCyclotomicPolynomial<Vector>(p::pq, p::Q), p::Q);

    // auto ct = TensorCt(ct0, ct1);

    // auto skk = TensorKey({sc0.skp}, {sc1.skp});

    // std::cout << "decrypted: " << sc0.RLWEDecrypt(ct, skk, p::t) << std::endl;

    // Scheme sctensor({p::ppq, pq});
    // Poly oneq = Poly(p::pp1, EVALUATION, true);
    // for (uint32_t i = 0; i < p::p1 - 1; i++) {
    //     oneq[i] = 1;
    // }
    // auto ksk = sctensor.KeySwitchGen(skk, {Tensor(sc0.skp, oneq)});
    // auto ctfinal = sctensor.KeySwitch(ct, ksk);

    return 0;
}