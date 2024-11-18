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
const usint N = 1024;
const usint Ncyc = 2 * N;
const usint p = 12289;
// const Integer Q("576460751449890817");
// const Integer g("5");
// const usint p = 769;
const usint rho = 4;
const usint Noverr = N / rho;

const Integer Q = lbcrypto::LastPrime<Integer>(60, p * (p - 1));
const Integer g = lbcrypto::FindGenerator(Q);
const Integer rou = g.ModExp((Q - 1) / p / (p - 1), Q);
const Integer rN = lbcrypto::RootOfUnity(Ncyc, Integer(p));
const Integer zeta = rN.ModExp(rho, p);
const usint Bks = 64;
const usint Rx = 8;

const auto nparams = std::make_shared<ILParams>(Ncyc, p, rN);
const auto pparams = std::make_shared<ILParams>(p, Q, rou, 0, 0);

}

Poly Monomial(usint k);

Vector PartialFourierTransform(Vector a, usint rho);
Vector PartialInverseFourierTransform(Vector a, usint rho, usint r);

RLWEGadgetCiphertext EvalInnerProduct(Scheme &sc, RLWEGadgetCiphertext ct, std::vector<RGSWCiphertext>::iterator z, const std::vector<Integer> &a, usint l, usint d = 1);
std::vector<RLWEGadgetCiphertext> HomomorphicPFT(Scheme &sc, std::vector<RLWEGadgetCiphertext> z);

int main() {
    primecyc::RaderFFTNat<Vector>::m_enabled[par::p] = true;

    BaseScheme base_sc({par::nparams, par::Ncyc, par::p, 2}, true);
    Scheme bt_sc({par::pparams, par::p, par::Q, par::Bks});

    base_sc.skp.SetFormat(COEFFICIENT);
    std::cout << base_sc.skp << std::endl;
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

    packct_a.SetFormat(COEFFICIENT);
    auto a_pft = PartialFourierTransform(packct_a.GetValues(), par::rho) * Integer(par::N / par::rho).ModInverse(par::p);
    packct_a.SetFormat(EVALUATION);

    Poly zero(par::pparams, EVALUATION, true);
    Poly one(par::pparams, EVALUATION, true);
    for (usint i = 0; i < par::p - 1; i++) {
        one[i] = 1;
    }
    RLWEGadgetCiphertext empty_reg;
    for (const auto &t : bt_sc.params.g) {
        empty_reg.push_back({zero, one * t});
    }

    std::vector<RLWEGadgetCiphertext> regs(par::N);
    Integer zzeta = par::zeta;
    for (usint k = 0, kk = 0; k < par::Noverr; k++) {
        for (usint i = 0; i < par::rho; i++) {
            std::vector<Integer> a(par::rho);
            for (usint j = 0; j <= i; j++) {
                a[i - j] = a_pft[kk * par::rho + j];
            }
            for (usint j = i + 1; j < par::rho; j++) {
                a[par::rho - j + i] = a_pft[kk * par::rho + j].ModMul(zzeta, par::p);
            }
            regs[kk * par::rho + i] = EvalInnerProduct(bt_sc, empty_reg, bk.begin() + kk * par::rho, a, par::rho);
        }
        zzeta.ModMulEq(par::zeta, par::p);
        zzeta.ModMulEq(par::zeta, par::p);
        for (usint l = par::Noverr >> 1; l > (kk ^= l); l >>= 1);
    }

    Poly foo = packct_a * base_sc.skp;
    foo.SetFormat(COEFFICIENT);
    std::cout << foo << std::endl;

    auto result = HomomorphicPFT(bt_sc, regs);
    for (usint i = 0; i < par::N; i++) {
        auto t = bt_sc.RLWEDecrypt(result[i], {bt_sc.skp}, par::t);
        usint tt = 0;
        if (t[0] != 0 && t[1] != 0) {
            tt = par::p - 1;
        } else {
            for (usint j = 0; j < par::p - 1; j++) {
                if (t[j] != 0) {
                    tt = j;
                    break;
                }
            }
        }
        std::cout << tt << " ";
    }
    std::cout << std::endl;

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

RLWEGadgetCiphertext EvalInnerProduct(Scheme &sc, RLWEGadgetCiphertext ct, std::vector<RGSWCiphertext>::iterator z, const std::vector<Integer> &a, usint l, usint d) {
    Integer t = 1;
    for (usint i = 0; i < l; i++) {
        if (a[i] != 0) {
            t.ModMulEq(a[i].ModInverse(par::p), par::p);
            if (t != 1) {
                ct = sc.GaloisConjugate(ct, t.ConvertToInt());
                ct = sc.KeySwitch(ct, sc.ksk_galois[t.ConvertToInt()]);
            }
            ct = sc.ExtMult(ct, z[i * d]);
            t = a[i];
        }
    }
    if (t != 1) {
        ct = sc.GaloisConjugate(ct, t.ConvertToInt());
        ct = sc.KeySwitch(ct, sc.ksk_galois[t.ConvertToInt()]);
    }
    return ct;
}

std::vector<RLWEGadgetCiphertext> HomomorphicPFT(Scheme &sc, std::vector<RLWEGadgetCiphertext> z) {
    std::vector<RGSWCiphertext> regs(par::N);
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
            Integer Z = 1;
            for (usint l = 0; l < i; l += par::rho) {
                for (usint m = 0; m < par::rho; m++) {
                    Integer zz = Z;
                    for (usint f = 0; f < r; f++) {
                        Integer zzz = zz;
                        std::vector<Integer> a(r - 1);
                        for (usint g = 1; g < r; g++) {
                            a[index[g] - 1] = zzz;
                            zzz.ModMulEq(zz, par::p);
                        }
                        z[k + l + m + f * i] = EvalInnerProduct(sc, regs[k + l + m].second, regs.begin() + k + l + m + i, a, r - 1, i);
                        zz.ModMulEq(w, par::p);
                    }
                }
                Z.ModMulEq(W, par::p);
            }
        }
    }

    Integer w = par::rN.ModExp(par::Ncyc - par::rho, par::p);
    Integer Z = 1;
    for (usint i = 0; i < par::N; i += par::rho) {
        for (usint j = 0; j < par::rho; j++) {
            if (Z != 1) {
                z[i + j] = sc.GaloisConjugate(z[i + j], Z.ConvertToInt());
                z[i + j] = sc.KeySwitch(z[i + j], sc.ksk_galois[Z.ConvertToInt()]);
            }
        }
        Z.ModMulEq(w, par::p);
    }

    return std::move(z);
}