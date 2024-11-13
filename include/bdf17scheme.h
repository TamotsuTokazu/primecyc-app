#ifndef _BDF17SCHEME_H_
#define _BDF17SCHEME_H_

#include "rlwe.h"

template <typename Poly>
class BDF17SchemeImpl : public SchemeImpl<Poly> {
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

    std::vector<RGSWCiphertext> bk;
    Poly x1;

    BDF17SchemeImpl(Params p, Poly x1_, Vector sk = {});

    RLWECiphertext RLWEEncrypt(const Poly &m, const RLWEKey &sk, const Integer &q_plain) override;

    RLWECiphertext Process(Vector a, Integer b, Integer q_plain);
};

#include "bdf17scheme-impl.h"

#endif // _BDF17SCHEME_H_
