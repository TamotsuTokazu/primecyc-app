#include "binfhecontext.h"
#include "params.h"
#include "rlwe.h"
#include "math/math-hal.h"

int main() {
    std::cout << p::nttSizepq << std::endl;

    uint32_t pq = p::p0 * p::p1;

    primecyc::RaderFFTNat<Vector>::m_enabled[p::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[p::p1] = true;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(pq);
    Vector sk = dugpq.GenerateVector(p::n);
    Vector a = dugpq.GenerateVector(p::n);
    for (uint32_t xx = 0; xx < p::pq; xx++) {

        Integer b = xx;
        for (uint32_t i = 0; i < p::n; i++) {
            b += a[i] * sk[i];
        }

        Scheme sc0{{p::pp0, p::p0}, sk};
        auto ct0 = sc0.Process(a, b, p::t);
        // std::cout << "ct0: " << sc0.RLWEDecrypt(ct0, {sc0.skp}, p::t) << std::endl;

        Scheme sc1{{p::pp1, p::p1}, sk};
        auto ct1 = sc1.Process(a, b, p::t);
        // std::cout << "ct1: " << sc1.RLWEDecrypt(ct1, {sc1.skp}, p::t) << std::endl;

        ChineseRemainderTransformArb<Vector>().SetCylotomicPolynomial(lbcrypto::GetCyclotomicPolynomial<Vector>(p::pq, p::Q), p::Q);

        auto ct = TensorCt(ct0, ct1);
        auto skk = TensorKey({sc0.skp}, {sc1.skp});

        // std::cout << "decrypted: " << sc0.RLWEDecrypt(ct, skk, p::t) << std::endl;

        Scheme scheme_tensor({p::ppq, pq});
        Poly skpq = Poly(p::ppq, COEFFICIENT, true);
        for (uint32_t i = 0; i < p::n; i++) {
            skpq[i * p::p1] = sk[i];
        }
        // std::cout << "skpq: " << skpq << std::endl;
        skpq.SetFormat(EVALUATION);
        auto tensor_ksk = scheme_tensor.KeySwitchGen(skk, {skpq});
        auto tensor_ct = scheme_tensor.KeySwitch(ct, tensor_ksk);

        auto plain = scheme_tensor.RLWEDecrypt(tensor_ct, {skpq}, p::Q);
        // std::cout << "decrypted: " << plain << std::endl;
        auto temp = TracePqToP(plain);
        // std::cout << "decrypted: " << temp << std::endl;
        scheme_tensor.ModSwitch(plain, p::t);
        // std::cout << "decrypted: " << plain << std::endl;
        scheme_tensor.ModSwitch(temp, p::t);
        // std::cout << "decrypted: " << temp << std::endl;

        skpq.SetFormat(COEFFICIENT);
        auto cta = tensor_ct[0];
        auto ctb = tensor_ct[1];
        cta.SetFormat(COEFFICIENT);
        ctb.SetFormat(COEFFICIENT);
        cta = TracePqToP(cta);
        ctb = TracePqToP(ctb);
        skpq = TracePqToP(skpq) * Integer(p::p1 - 1).ModInverse(p::Q);
        // std::cout << "cta: " << cta << std::endl;
        // std::cout << "skpq: " << skpq << std::endl;
        cta.SetFormat(EVALUATION);
        ctb.SetFormat(EVALUATION);
        skpq.SetFormat(EVALUATION);

        auto plain1 = scheme_tensor.RLWEDecrypt({cta, ctb}, {skpq}, p::Q);
        // std::cout << "decrypted: " << plain1 << std::endl;
        // scheme_tensor.ModSwitch(plain1, p::t);
        // std::cout << "decrypted: " << plain1 << std::endl;

        Integer z = TracePtoZ(plain1);
        // std::cout << "decrypted: " << z << std::endl;
        z = (z * p::t + (p::Q / 2)) / p::Q;
        std::cout << "decrypted: " << z << std::endl;

    }

    return 0;
}