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
    sc0.Process(a, b);

    Scheme sc1{{p::pp1, p::p1}, sk};
    sc1.Process(a, b);

    return 0;
}