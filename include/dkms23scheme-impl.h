#ifndef _DKMS23SCHEME_IMPL_H_
#define _DKMS23SCHEME_IMPL_H_

template <typename Poly>
DKMS23SchemeImpl<Poly>::DKMS23SchemeImpl(Params p) : SchemeImpl<Poly>(p, true) {
    this->GaloisKeyGen();
}

#endif // _DKMS23SCHEME_IMPL_H_