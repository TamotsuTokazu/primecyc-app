#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

using Integer = lbcrypto::NativeInteger;
using Vector = lbcrypto::NativeVector;
using Poly = lbcrypto::NativePoly;

int main() {
    usint p = 5;
    usint q = 7;
    usint pq = p * q;
    primecyc::TensorFFTNat<Vector>::m_factor[pq] = {p, q};
    Integer Q = lbcrypto::FirstPrime<Integer>(56, pq * (p - 1) * (q - 1));
    Integer r = lbcrypto::RootOfUnity<Integer>(pq * (p - 1) * (q - 1), Q);
    Integer rp = r.ModExp(q * (q - 1), Q);
    Integer rq = r.ModExp(p * (p - 1), Q);
    // Integer Q = 36066736134770689;
    // Integer r = 4364918564594134;
    usint tot = (p - 1) * (q - 1);

    Vector ap(p - 1, Q);
    Vector aq(q - 1, Q);

    // ap[3] = 1;
    for (usint i = 0; i < p - 1; i++) {
        ap[i] = Q - 1;
    }
    for (usint j = 0; j < q - 1; j++) {
        aq[j] = Q - 1;
    }
    // aq[4] = 1;

    Vector a(tot, Q);

    // ap[0] = Integer(0).ModSub(ap[0], Q);
    // for (usint i = 1; i < p - 1; i++) {
    //     ap[i].ModAddEq(ap[0], Q);
    // }

    // aq[0] = Integer(0).ModSub(aq[0], Q);
    // for (usint i = 1; i < q - 1; i++) {
    //     aq[i].ModAddEq(aq[0], Q);
    // }

    // std::cout << "ap: " << ap << std::endl;
    // std::cout << "aq: " << aq << std::endl;

    ap = primecyc::RaderFFTNat<Vector>().ForwardRaderPermute(ap, rp);
    aq = primecyc::RaderFFTNat<Vector>().ForwardRaderPermute(aq, rq);

    // std::cout << "ap: " << ap << std::endl;
    // std::cout << "aq: " << aq << std::endl;

    for (usint i = 0; i < p - 1; i++) {
        for (usint j = 0; j < q - 1; j++) {
            a[(q - 1) * i + j] = aq[j].ModMul(ap[i], Q);
        }
    }

    a = primecyc::TensorFFTNat<Vector>().InverseTransform(a, r, pq);

    std::cout << "a: " << a << std::endl;

    return 0;
}