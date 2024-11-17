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
const usint N = 16;
const usint Ncyc = 2 * N;
const usint p = 12289;
// const Integer Q("576460751449890817");
// const Integer g("5");
// const usint p = 73;
const usint rho = 4;

const Integer Q = lbcrypto::LastPrime<Integer>(60, p * (p - 1));
const Integer g = lbcrypto::FindGenerator(Q);
const Integer rou = g.ModExp((Q - 1) / p / (p - 1), Q);
const Integer rN = lbcrypto::RootOfUnity(Ncyc, Integer(p));
const Integer zeta = rN.ModExp(2 * rho, p);
const usint Bks = 2;
const usint Rx = 2;

const auto nparams = std::make_shared<ILParams>(Ncyc, p, rN);
const auto pparams = std::make_shared<ILParams>(p, Q, rou, 0, 0);

}

Poly Monomial(usint k);

Vector PartialFourierTransform(Vector a, usint rho);
Vector PartialInverseFourierTransform(Vector a, usint rho, usint r);

RLWEGadgetCiphertext EvalInnerProduct(Scheme &sc, std::vector<RGSWCiphertext>::iterator z, const std::vector<Integer> &a);
std::vector<RLWEGadgetCiphertext> HomomorphicPFT(Scheme &sc, std::vector<RLWEGadgetCiphertext> z);

int main() {
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p] = true;

    Vector a(par::N, par::p);
    for (usint i = 0; i < par::N; i++) {
        a[i] = 1;
        std::cout << a << std::endl;
        Vector b = PartialFourierTransform(a, par::rho);
        std::cout << b << std::endl;
        Vector c = PartialInverseFourierTransform(b, par::rho, par::Rx);
        std::cout << c << std::endl << std::endl;
        a[i] = 0;
    }
    return 0;

    BaseScheme base_sc({par::nparams, par::Ncyc, par::p, 2}, true);
    Scheme bt_sc({par::pparams, par::p, par::Q, par::Bks});

    Vector z_pft = PartialFourierTransform(base_sc.skp.GetValues(), par::rho);
    std::cout << z_pft << std::endl;
    std::vector<RGSWCiphertext> bk;
    for (usint i = 0; i < par::N; i++) {
        bk.push_back(bt_sc.RGSWEncrypt(Monomial(z_pft[i].ConvertToInt()), {bt_sc.skp}));
    }

    base_sc.skp.SetFormat(EVALUATION);

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
    Poly packct_b(par::nparams, COEFFICIENT, true);
    for (usint i = 0; i < par::N; i++) {
        packct_b[i] = lwe_b[i];
    }
    packct_b.SetFormat(EVALUATION);
    for (usint i = 0; i < par::n; i++) {
        Poly t(par::nparams, COEFFICIENT, true);
        for (usint j = 0; j < par::N; j++) {
            t[j] = lwe_a[j][i];
        }
        auto ct = base_sc.Mult(t, packing_key[i]);
        packct_a -= ct[0];
        packct_b -= ct[1];
    }
    RLWECiphertext pack_ct = {packct_a, packct_b};

    auto a_pft = PartialFourierTransform(packct_a.GetValues(), par::rho);

    Poly zero(par::pparams, EVALUATION, true);
    Poly one(par::pparams, EVALUATION, true);
    for (usint i = 0; i < par::p - 1; i++) {
        zero[i] = 1;
    }

    std::vector<RLWEGadgetCiphertext> regs(par::N);
    for (usint l = 0; l < par::N; l += par::rho) {
        for (usint i = 0; i < par::rho; i++) {
            std::vector<Integer> a(par::rho);
            for (usint j = 0; j <= i; j++) {
                a[i - j] = a_pft[j];
                std::cout << i - j << std::endl;
            }
            for (usint j = i + 1; j < par::rho; j++) {
                a[par::rho - j + i] = a_pft[j].ModMul(par::zeta, par::Q);
                std::cout << par::rho - j + i << std::endl;
            }
            regs[l + i] = EvalInnerProduct(bt_sc, bk.begin() + l, a);
        }
    }

    std::vector<RLWEGadgetCiphertext> regs2(par::N);

    for (usint i = 1; par::rho * i <= par::N; i *= par::Rx) {
        for (usint j = 0; j < par::N / i; j += par::rho) {
            std::vector<RGSWCiphertext> z;
            for (usint k = 0; k < par::rho; k++) {
                z.push_back(bt_sc.RGSWEncrypt(Monomial(j + k), {bt_sc.skp}));
            }
            for (usint k = 0; k < par::rho; k++) {
                regs2[j + k] = bt_sc.ExtMult(regs2[j + k], z[k]);
            }
        }
    }
    return 0;
}

Poly Monomial(usint k) {
    Poly result(par::pparams, COEFFICIENT, true);
    if (k == par::p - 1) {
        for (usint i = 0; i < par::p - 1; i++) {
            result[i] = par::Q - 1;
        }
    } else {
        result[k] = 1;
    }
    result.SetFormat(EVALUATION);
    return result;
}

Vector PartialFourierTransform(Vector a, usint rho) {
    usint n = a.GetLength();
    for (usint i = n; i > rho; i >>= 1) {
        usint j = i >> 1;
        usint nn = n / i;
        Integer z = par::rN.ModExp(j, par::p);
        Integer w = par::rN.ModExp(i, par::p);
        for (usint k = 0, kk = 0; k < nn; k++) {
            for (usint l = 0; l < j; l++) {
                Integer u = a[kk * i + l + j].ModMul(z, par::p);
                a[kk * i + l + j] = a[kk * i + l].ModSub(u, par::p);
                a[kk * i + l].ModAddEq(u, par::p);
            }
            for (usint l = nn >> 1; l > (kk ^= l); l >>= 1);
            z.ModMulEq(w, par::p);
        }
    }
    return std::move(a);
}

Vector PartialInverseFourierTransform(Vector a, usint rho, usint r) {
    usint n = a.GetLength();
    Vector temp(r);
    std::vector<usint> index;

    for (usint i = rho; i < n; i *= r) {
        usint j = std::min(i * r, n);
        r = j / i;
        Integer W = par::rN.ModExp(par::Ncyc - par::Ncyc / (j / rho), par::p);
        Integer w = par::rN.ModExp(par::Ncyc - par::Ncyc / r, par::p);
        if (r != index.size()) {
            index.resize(r);
            for (usint k = 0, kk = 0; k < r; k++) {
                index[k] = kk;
                for (usint l = r >> 1; l > (kk ^= l); l >>= 1);
            }
        }
        for (usint k = 0; k < n; k += j) {
            Integer z = 1;
            for (usint l = 0; l < i; l += rho) {
                for (usint m = 0; m < rho; m++) {
                // first group: a[k + m], a[k + m + i], ..., a[k + m + (r - 1) * i]
                // second group: a[k + m + rho], a[k + m + rho + i], ..., a[k + m + rho + (r - 1) * i]
                // number of groups: i / rho
                // number of elements in each group: r
                    Integer zz = z;
                    for (usint f = 0; f < r; f++) {
                        auto &t = temp[f];
                        Integer zzz = zz;
                        t = a[k + l + m];
                        for (usint g = 1; g < r; g++) {
                            t.ModAddEq(a[k + l + m + index[g] * i].ModMul(zzz, par::p), par::p);
                            zzz.ModMulEq(zz, par::p);
                        }
                        zz.ModMulEq(w, par::p);
                    }
                    for (usint f = 0; f < r; f++) {
                        a[k + l + m + f * i] = temp[f];
                    }
                }
                z.ModMulEq(W, par::p);
            }
        }
    }
    Integer w = par::rN.ModExp(par::Ncyc - rho, par::p);
    Integer z = 1;
    for (usint i = 0; i < n; i += rho) {
        for (usint j = 0; j < rho; j++) {
            a[i + j].ModMulEq(z, par::p);
        }
        z.ModMulEq(w, par::p);
    }
    return std::move(a);
}

RLWEGadgetCiphertext EvalInnerProduct(Scheme &sc, RLWEGadgetCiphertext ct, std::vector<RGSWCiphertext>::iterator z, const std::vector<Integer> &a) {
    Integer t = 1;
    for (usint i = 0; i < par::rho; i++) {
        if (a[i] != 0) {
            t.ModMulEq(a[i].ModInverse(par::p), par::p);
            if (t != 1) {
                ct = sc.GaloisConjugate(ct, t.ConvertToInt());
                ct = sc.KeySwitch(ct, sc.ksk_galois[t.ConvertToInt()]);
            }
            ct = sc.ExtMult(ct, z[i]);
        }
    }
    if (t != 1) {
        ct = sc.GaloisConjugate(ct, t.ConvertToInt());
        ct = sc.KeySwitch(ct, sc.ksk_galois[t.ConvertToInt()]);
    }
    std::cout << "decrypted: " << sc.RLWEDecrypt(ct, {sc.skp}, par::t) << std::endl;
    return ct;
}

std::vector<RLWEGadgetCiphertext> HomomorphicPFT(Scheme &sc, std::vector<RLWEGadgetCiphertext> z) {
    std::vector<RGSWCiphertext> regs;
    std::vector<usint> index;

    for (usint i = par::rho; i < par::N; i *= par::Rx) {
        usint j = std::min(i * par::Rx, par::N);
        usint r = j / i;
        Integer W = par::rN.ModExp(par::Ncyc - par::Ncyc / (j / par::rho), par::p);
        Integer w = par::rN.ModExp(par::Ncyc - par::Ncyc / r, par::p);
        if (r != index.size()) {
            index.resize(r);
            for (usint k = 0, kk = 0; k < r; k++) {
                index[k] = kk;
                for (usint l = r >> 1; l > (kk ^= l); l >>= 1);
            }
        }
        for (usint k = 0; k < par::N; k++) {
            regs[k] = sc.SchemeSwitch(z[k]);
        }
        for (usint k = 0; k < par::N; k += j) {
            Integer z = 1;
            for (usint l = 0; l < i; l += par::rho) {
                for (usint m = 0; m < rho; m++) {
                // first group: a[k + m], a[k + m + i], ..., a[k + m + (r - 1) * i]
                // second group: a[k + m + rho], a[k + m + rho + i], ..., a[k + m + rho + (r - 1) * i]
                // number of groups: i / rho
                // number of elements in each group: r
                    Integer zz = z;
                    for (usint f = 0; f < r; f++) {
                        auto &t = temp[f];
                        Integer zzz = zz;
                        t = a[k + l + m];
                        for (usint g = 1; g < r; g++) {
                            t.ModAddEq(a[k + l + m + index[g] * i].ModMul(zzz, par::p), par::p);
                            zzz.ModMulEq(zz, par::p);
                        }
                        zz.ModMulEq(w, par::p);
                    }
                    for (usint f = 0; f < r; f++) {
                        a[k + l + m + f * i] = temp[f];
                    }
                }
                z.ModMulEq(W, par::p);
            }
        }
    }
    Integer w = par::rN.ModExp(par::Ncyc - rho, par::p);
    Integer z = 1;
    for (usint i = 0; i < n; i += rho) {
        for (usint j = 0; j < rho; j++) {
            a[i + j].ModMulEq(z, par::p);
        }
        z.ModMulEq(w, par::p);
    }

    return std::move(z);
}