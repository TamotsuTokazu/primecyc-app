#include "binfhecontext.h"
#include "dkms23scheme.h"
#include "math/math-hal.h"

#include <chrono>

using Poly = lbcrypto::NativePoly;
using Scheme = DKMS23SchemeImpl<Poly>;
using BaseScheme = SchemeImpl<Poly>;
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
const usint N = 4;
const usint Ncyc = 2 * N;
const usint p = 12289;
// const usint p = 13;
// const Integer Q = lbcrypto::LastPrime<Integer>(59, p * (p - 1));
const Integer Q("576460751449890817");
const Integer g("5");
const Integer rou = g.ModExp((Q - 1) / p / (p - 1), Q);
const Integer rN = lbcrypto::RootOfUnity(Ncyc, Integer(p));
const usint Bks = 256;

const auto nparams = std::make_shared<ILParams>(Ncyc, p, rN);
const auto pparams = std::make_shared<ILParams>(p, Q, rou, 0, 0);

}

int main() {
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p] = true;

    BaseScheme base_sc({par::nparams, par::Ncyc, par::p, 2}, true);
    // Scheme bt_sc({par::pparams, par::p, par::Q, par::Bks});

    Vector sk = lbcrypto::DiscreteGaussianGeneratorImpl<Vector>().GenerateVector(par::n, par::p);

    std::vector<RLWEGadgetCiphertext> packing_key;
    for (usint i = 0; i < par::n; i++) {
        Poly m(par::nparams, EVALUATION, true);
        for (usint j = 0; j < par::N; j++) {
            m[j] = sk[i];
        }
        packing_key.push_back(base_sc.RLWEGadgetEncrypt(m, {base_sc.skp}, par::p));
    }

    std::vector<Vector> lwe_a(par::N);
    std::vector<Integer> lwe_b(par::N);
    auto dug = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(par::p);
    for (usint i = 0; i < par::N; i++) {
        lwe_a[i] = dug.GenerateVector(par::n);
        lwe_b[i] = 0;
        for (usint j = 0; j < par::n; j++) {
            lwe_b[i].ModAddEq(lwe_a[i][j] * sk[j], par::p);
        }
    }

    Poly packct_a(par::nparams, EVALUATION, true);
    Poly packct_b(par::nparams, EVALUATION, true);
    for (usint i = 0; i < par::n; i++) {
        Poly t(par::nparams, COEFFICIENT, true);
        for (usint j = 0; j < par::N; j++) {
            t[j] = lwe_a[j][i];
        }
        auto ct = base_sc.Mult(t, packing_key[i]);
        packct_a -= ct[0];
        packct_b += ct[1];
    }
    RLWECiphertext pack_ct = {packct_a, packct_b};
    return 0;
}