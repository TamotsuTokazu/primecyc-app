#include "binfhecontext.h"
#include "math/math-hal.h"

using Vector = lbcrypto::NativeVector;
using Integer = lbcrypto::NativeInteger;

const usint N = 7;
const Integer Q = lbcrypto::LastPrime<Integer>(59, N * (N - 1));
const Integer rootOfUnity = lbcrypto::RootOfUnity(N * (N - 1), Q);

int main() {
    primecyc::RaderFFTNat<Vector>::m_enabled[N] = true;
    Vector a(N - 1, Q);
    a[0] = 1;
    Vector b = primecyc::RaderFFTNat<Vector>().ForwardRaderPermute(a, rootOfUnity);
    std::cout << b << std::endl;
    b = primecyc::RaderFFTNat<Vector>().InverseRaderPermute(b, rootOfUnity);
    std::cout << b << std::endl;
    return 0;
}