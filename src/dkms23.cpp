#include "binfhecontext.h"
#include "rlwe.h"
#include "math/math-hal.h"

#include <chrono>

using Poly = lbcrypto::NativePoly;
using Scheme = SchemeImpl<Poly>;
using Vector = Scheme::Vector;
using Integer = Scheme::Integer;
using ILParams = Scheme::ILParams;
using Params = Scheme::Params;

using RLWECiphertext = Scheme::RLWECiphertext;
using RLWEKey = Scheme::RLWEKey;
using RLWESwitchingKey = Scheme::RLWESwitchingKey;
using RGSWCiphertext = Scheme::RGSWCiphertext;

using lbcrypto::BigInteger;

namespace par {

const usint N = 1024;
const usint p = 12289;
// const Integer Q = lbcrypto::FirstPrime<Integer>(60, p * (p - 1));
const Integer Q("1152921505164890113");
const Integer g("5");
const Integer rOU = g.ModExp((Q - 1) / p / (p - 1), Q);

}

int main() {
    return 0;
}