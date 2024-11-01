#include "binfhecontext.h"
#include "params.h"
#include "rlwe-impl.h"

int main() {
    uint32_t pq = p::p0 * p::p1;
    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(pq);
    Vector sk = dugpq.GenerateVector(p::n);
    Vector a = dugpq.GenerateVector(p::n);
    Integer b = 1;
    for (uint32_t i = 0; i < p::n; i++) {
        b += a[i] * sk[i];
    }

    Scheme sc0{{p::pp0, p::p0}, sk};
    auto ct0 = sc0.Process(a, b, p::t);
    std::cout << "ct0: " << sc0.RLWEDecrypt(ct0, {sc0.skp}, p::t) << std::endl;

    Scheme sc1{{p::pp1, p::p1}, sk};
    auto ct1 = sc1.Process(a, b, p::t);
    std::cout << "ct1: " << sc1.RLWEDecrypt(ct1, {sc1.skp}, p::t) << std::endl;

    ChineseRemainderTransformArb<Vector>().SetCylotomicPolynomial(lbcrypto::GetCyclotomicPolynomial<Vector>(p::pq, p::Q), p::Q);

    auto ct = TensorCt(ct0, ct1);
    std::cout << "ct: " << ct << std::endl;

    auto skk = TensorKey({sc0.skp}, {sc1.skp});
    std::cout << "skk: " << skk << std::endl;

    std::cout << "decrypted: " << sc0.RLWEDecrypt(ct, skk, p::t) << std::endl;

    return 0;
}