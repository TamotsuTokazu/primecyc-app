#include "binfhecontext.h"
#include "math/math-hal.h"

#include <iostream>
#include <chrono>

#include <valgrind/callgrind.h>

using Poly = lbcrypto::NativePoly;
using Integer = Poly::Integer;
using ILParams = Poly::Params;

const Integer t("64");
const usint n = 60;
const usint p0 = 1153;
const usint p1 = 1297;
const usint pq = p0 * p1;
const usint Bks = 1 << 4;
const Integer Q("1152920977604149249");

const Integer rootOfUnity("739447795444923848");
const Integer rootOfUnity0("587824393635068054");
const Integer rootOfUnity1("950015890856494149");

const auto pp0 = std::make_shared<ILParams>(p0, Q, rootOfUnity0, 0, 0);

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    primecyc::RaderFFTNat<Poly::Vector>::m_enabled[p0] = true;
    CALLGRIND_START_INSTRUMENTATION;
    for (int i = 0; i < 100000; i++) {
        Poly a(pp0, COEFFICIENT, true);
        Poly b(pp0, COEFFICIENT, true);
        a[0] = i + 1;
        b[0] = i;
        a.SetFormat(EVALUATION);
        b.SetFormat(EVALUATION);
        Poly c = a * b;
        c.SetFormat(COEFFICIENT);
    }
    CALLGRIND_STOP_INSTRUMENTATION;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    return 0;
}