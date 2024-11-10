#include "binfhecontext.h"
#include "math/math-hal.h"

#include <chrono>

using Vector = lbcrypto::NativeVector;
using Integer = Vector::Integer;

const usint p = 1153;
const Integer Q("576485048298018817");

const Integer w("397835389725933388");

int main() {
    auto start = std::chrono::system_clock::now();

    Vector a(p - 1, Q);
    Vector b(p - 1, Q);

    const usint N = 100000;
    primecyc::RaderFFTNat<Vector>::m_enabled[p] = true;

    for (usint i = 0; i < N; i++) {
        a[0] = i;
        b[0] = i + 1;
        Vector at = primecyc::RaderFFTNat<Vector>().ForwardRaderPermute(a, w);
        Vector bt = primecyc::RaderFFTNat<Vector>().ForwardRaderPermute(b, w);
        Vector c = at * bt;
        Vector d = primecyc::RaderFFTNat<Vector>().InverseRaderPermute(c, w);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Time: " << elapsed_seconds.count() << std::endl;

    return 0;
}