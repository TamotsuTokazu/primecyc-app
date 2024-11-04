#include "binfhecontext.h"
#include "params.h"
#include "rlwe.h"
#include "math/math-hal.h"

int main() {
    uint32_t pq = p::p0 * p::p1;

    primecyc::RaderFFTNat<Vector>::m_enabled[p::p0] = true;
    primecyc::RaderFFTNat<Vector>::m_enabled[p::p1] = true;

    auto dugpq = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(pq);
    Vector sk = dugpq.GenerateVector(p::n);
    Vector a = dugpq.GenerateVector(p::n);

    Integer b = 2000;
    for (uint32_t i = 0; i < p::n; i++) {
        b += a[i] * sk[i];
    }

    // measure the running time the following code
    Scheme sc0{{p::pp0, p::p0}, sk};
    Scheme sc1{{p::pp1, p::p1}, sk};

    auto start = std::chrono::high_resolution_clock::now();

    auto ct0 = sc0.Process(a, b, p::t);
    std::cout << "ct0: " << sc0.RLWEDecrypt(ct0, {sc0.skp}, p::t) << std::endl;

    auto ct1 = sc1.Process(a, b, p::t);
    std::cout << "ct1: " << sc1.RLWEDecrypt(ct1, {sc1.skp}, p::t) << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;


    return 0;
}