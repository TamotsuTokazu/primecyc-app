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

    DKMS23SchemeImpl(Params p);
};

#include "dkms23scheme-impl.h"

#endif // _DKMS23SCHEME_H_