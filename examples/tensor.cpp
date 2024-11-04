#include "binfhecontext.h"
#include "params.h"
#include "rlwe.h"
#include "math/math-hal.h"

using Integer = lbcrypto::NativeInteger;
using Vector = lbcrypto::NativeVector;
using Poly = lbcrypto::NativePoly;

int main() {
    usint p = 1153;
    usint q = 1297;
    usint pq = p * q;
    primecyc::TensorFFTNat<Vector>::m_factor[pq] = {p, q};
    // Integer Q = lbcrypto::FirstPrime<Integer>(56, pq * (p - 1) * (q - 1));
    // Integer r = lbcrypto::RootOfUnity<Integer>(pq * (p - 1) * (q - 1), Q);
    Integer Q = 36066736134770689;
    Integer r = 4364918564594134;
    usint tot = (p - 1) * (q - 1);
    Vector a(tot, Q);
    a[0] = 1;
    // std::cout << "a: " << a << std::endl;
    Vector b = primecyc::TensorFFTNat<Vector>().ForwardTransform(a, r, pq);
    std::cout << "b: " << b << std::endl;

    Vector c = primecyc::TensorFFTNat<Vector>().InverseTransform(b, r, pq);
    // std::cout << "c: " << c << std::endl;
    return 0;
}