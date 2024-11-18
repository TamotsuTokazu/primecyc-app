#ifndef _DKMS23SCHEME_H_
#define _DKMS23SCHEME_H_

#include "rlwe.h"

template <typename Poly>
class DKMS23SchemeImpl : public SchemeImpl<Poly> {
public:
    using Vector = typename SchemeImpl<Poly>::Vector;
    using Integer = typename SchemeImpl<Poly>::Integer;
    using ILParams = typename SchemeImpl<Poly>::ILParams;

    using RLWEKey = typename SchemeImpl<Poly>::RLWEKey;
    using RLWECiphertext = typename SchemeImpl<Poly>::RLWECiphertext;
    using RLWEGadgetCiphertext = typename SchemeImpl<Poly>::RLWEGadgetCiphertext;
    using RLWESwitchingKey = typename SchemeImpl<Poly>::RLWESwitchingKey;
    using RGSWCiphertext = typename SchemeImpl<Poly>::RGSWCiphertext;

    using Params = typename SchemeImpl<Poly>::Params;

    RLWEGadgetCiphertext ssk;

    DKMS23SchemeImpl(Params p);

    Poly RLWEDecrypt(const RLWEGadgetCiphertext &ct, const RLWEKey &sk, const Integer &q_plain);
    RLWEGadgetCiphertext KeySwitch(const RLWEGadgetCiphertext &ct, const RLWESwitchingKey &ksk);
    RLWEGadgetCiphertext ExtMult(const RLWEGadgetCiphertext &ct, const RGSWCiphertext &ctGSW);
    RGSWCiphertext SchemeSwitch(const RLWEGadgetCiphertext &ct);
};

#include "dkms23scheme-impl.h"

#endif // _DKMS23SCHEME_H_