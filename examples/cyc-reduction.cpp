#include "binfhecontext.h"
#include "math/math-hal.h"

using Vector = lbcrypto::BigVector;
using Integer = Vector::Integer;

int main() {
    usint p = 3;
    usint q = 5;
    usint pq = p * q;
    usint cycloOrder = pq;
    // usint n = (p - 1) * (q - 1);
    Integer modulus = 127;
    Vector element(pq, modulus);
    element[13] = 1;
    element[14] = 1;
    auto l = element.GetLength();
    Vector a(l + p + q, modulus);
    std::cout << "here " << element << std::endl;
    for (usint i = 0; i < l; i++) {
        a[i] = element[i];
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = a.GetLength() - 1; i >= p; i--) {
        a[i].ModAddEq(a[i - p], modulus);
        a[i - p] = Integer(0).ModSub(a[i - p], modulus);
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = a.GetLength() - 1; i >= q; i--) {
        a[i].ModAddEq(a[i - q], modulus);
        a[i - q] = Integer(0).ModSub(a[i - q], modulus);
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = 0; i < l + p + q - 1; i++) {
        a[i] = Integer(0).ModSub(a[i], modulus);
        a[i + 1].ModSubEq(a[i], modulus);
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = a.GetLength() - 1; i >= cycloOrder; i--) {
        a[i - cycloOrder].ModAddEq(a[i], modulus);
        a[i] = 0;
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = a.GetLength() - 1; i >= 1; i--) {
        a[i].ModAddEq(a[i - 1], modulus);
        a[i - 1] = Integer(0).ModSub(a[i - 1], modulus);
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = 0; i < l + p; i++) {
        a[i] = Integer(0).ModSub(a[i], modulus);
        a[i + q].ModSubEq(a[i], modulus);
    }
    std::cout << "a: " << a << std::endl;
    for (usint i = 0; i < l + q; i++) {
        a[i] = Integer(0).ModSub(a[i], modulus);
        a[i + p].ModSubEq(a[i], modulus);
    }
    std::cout << "a: " << a << std::endl;
    return 0;
}