#ifndef _RLWE_H_
#define _RLWE_H_

#include "binfhecontext.h"
#include <vector>
#include "params.h"

struct Scheme {

    Params params;

    std::vector<RLWESwitchingKey> ksk_galois;
    std::vector<RGSWCiphertext> bk;

    Poly skp;

    Scheme(Params p, Vector sk = {});

    Poly GaloisConjugate(const Poly &x, const uint32_t &a);

    template <typename T>
    std::vector<T> GaloisConjugate(const std::vector<T> &x, const uint32_t &a);

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain);
    Poly RLWEDecrypt(const RLWECiphertext &ct, const RLWEKey &sk, const Integer &q_plain);

    void ModSwitch(Poly &x, const Integer &q);
    void ModSwitch(RLWECiphertext &ct, const Integer &q);

    RLWESwitchingKey KeySwitchGen(const RLWEKey &sk, const RLWEKey &skN);

    RLWECiphertext KeySwitch(const RLWECiphertext &ct, const RLWESwitchingKey &K);

    RGSWCiphertext RGSWEncrypt(const Poly &m, const RLWEKey &sk);

    RLWECiphertext ExtMult(const RLWECiphertext &ct, const RGSWCiphertext &ctGSW);

    RLWECiphertext Process(Vector a, Integer b, Integer q_plain);
};

Poly Tensor(const Poly &a, const Poly &b);

RLWEKey TensorKey(const RLWEKey &skp, const RLWEKey &skq);

RLWECiphertext TensorCt(const RLWECiphertext &ap, const RLWECiphertext &aq);

#include "rlwe-impl.h"

#endif