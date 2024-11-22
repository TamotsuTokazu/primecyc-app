#ifndef _DKMS23SCHEME_IMPL_H_
#define _DKMS23SCHEME_IMPL_H_

template <typename Poly>
DKMS23SchemeImpl<Poly>::DKMS23SchemeImpl(Params p) : SchemeImpl<Poly>(p, true) {
    this->GaloisKeyGen();
    ssk = this->RLWEGadgetEncrypt(this->skp * this->skp, {this->skp}, this->params.Q);
}

template <typename Poly>
Poly DKMS23SchemeImpl<Poly>::RLWEDecrypt(const RLWEGadgetCiphertext &ct, const RLWEKey &sk, const Integer &q_plain) {
    Poly m(this->params.poly, COEFFICIENT, true);
    m[0] = this->params.Q / q_plain;
    return SchemeImpl<Poly>::RLWEDecrypt(this->Mult(m, ct), sk, q_plain);
}

template <typename Poly>
typename DKMS23SchemeImpl<Poly>::RLWEGadgetCiphertext DKMS23SchemeImpl<Poly>::KeySwitch(const RLWEGadgetCiphertext &ct, const RLWESwitchingKey &ksk) {
    RLWEGadgetCiphertext result;
    for (const auto &cti: ct) {
        result.push_back(SchemeImpl<Poly>::KeySwitch(cti, ksk));
    }
    return result;
}

template <typename Poly>
typename DKMS23SchemeImpl<Poly>::RLWEGadgetCiphertext DKMS23SchemeImpl<Poly>::ExtMult(const RLWEGadgetCiphertext &ct, const RGSWCiphertext &ctGSW) {
    RLWEGadgetCiphertext result;
    for (const auto &cti: ct) {
        result.push_back(SchemeImpl<Poly>::ExtMult(cti, ctGSW));
    }
    return result;
}

template <typename Poly>
typename DKMS23SchemeImpl<Poly>::RGSWCiphertext DKMS23SchemeImpl<Poly>::SchemeSwitch(const RLWEGadgetCiphertext &ct) {
    RLWEGadgetCiphertext ct2;
    for (const auto &cti: ct) {
        const auto &a = cti[0];
        const auto &b = cti[1];
        auto x = this->Mult(a.Negate(), ssk);
        x[0] -= b;
        ct2.push_back(x);
    }
    return {ct2, ct};
}

#endif // _DKMS23SCHEME_IMPL_H_