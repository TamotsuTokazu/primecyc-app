#include "binfhecontext.h"
#include "dkms23scheme.h"
#include "math/math-hal.h"

#include <chrono>

using Poly = lbcrypto::NativePoly;
using Scheme = DKMS23SchemeImpl<Poly>;
using Vector = Scheme::Vector;
using Integer = Scheme::Integer;
using ILParams = Scheme::ILParams;
using Params = Scheme::Params;

using RLWECiphertext = Scheme::RLWECiphertext;
using RLWEGadgetCiphertext = Scheme::RLWEGadgetCiphertext;
using RLWEKey = Scheme::RLWEKey;
using RLWESwitchingKey = Scheme::RLWESwitchingKey;
using RGSWCiphertext = Scheme::RGSWCiphertext;

using lbcrypto::BigInteger;

namespace par {

const Integer t("4");
const usint n = 5;
const usint N = 512;
const usint Ncyc = 1024;
const usint p = 12289;
// const Integer Q = lbcrypto::LastPrime<Integer>(59, p * (p - 1));
const Integer Q("576460751449890817");
const Integer g("5");
const Integer rOU = g.ModExp((Q - 1) / p / (p - 1), Q);
const usint Bks = 256;

const auto pparams = std::make_shared<ILParams>(p, Q, rOU, 0, 0);

}

int main() {
    std::cout << "Q: " << par::Q << std::endl;
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p] = true;

    Scheme scheme({par::pparams, par::p, par::Q, par::Bks});

    Vector sk = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(par::p).GenerateVector(par::n);

    std::vector<RLWEGadgetCiphertext> packing_key;
    for (usint i = 0; i < par::n; i++) {
        Poly m(par::pparams, EVALUATION, true);
        for (usint j = 0; j < par::N; j++) {
            m[j] = sk[i];
        }
        packing_key.push_back(scheme.RLWEGadgetEncrypt(m, {scheme.skp}, par::Q));
    }
    return 0;
}