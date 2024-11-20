#include <benchmark/benchmark.h>
#include "math/math-hal.h"
#include "binfhecontext.h"

static void BM_PowerOf2FFT512(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 512;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, 2 * N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);

    for (auto _ : state) {
        Vector b(N, Q);
        ChineseRemainderTransformFTT<Vector>().ForwardTransformToBitReverse(a, rootOfUnity, N * 2, &b);
        a = b;
    }
}

static void BM_PowerOf2FFT1024(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 1024;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, 2 * N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);

    for (auto _ : state) {
        Vector b(N, Q);
        ChineseRemainderTransformFTT<Vector>().ForwardTransformToBitReverse(a, rootOfUnity, N * 2, &b);
        a = b;
    }
}

static void BM_PowerOf2FFT8192(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 8192;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, 2 * N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);

    for (auto _ : state) {
        Vector b(N, Q);
        ChineseRemainderTransformFTT<Vector>().ForwardTransformToBitReverse(a, rootOfUnity, N * 2, &b);
        a = b;
    }
}

static void BM_PowerOf2FFT16384(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 16384;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, 2 * N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);

    for (auto _ : state) {
        Vector b(N, Q);
        ChineseRemainderTransformFTT<Vector>().ForwardTransformToBitReverse(a, rootOfUnity, N * 2, &b);
        a = b;
    }
}

static void BM_ArbFFT769(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[769] = true;

    const usint N = 769;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, N * (N - 1));
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N * (N - 1), Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_ArbFFT1153(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[1153] = true;

    const usint N = 1153;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, N * (N - 1));
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N * (N - 1), Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_ArbFFT3889(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[3889] = true;

    const usint N = 3889;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, N * (N - 1));
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N * (N - 1), Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_ArbFFT10369(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[10369] = true;

    const usint N = 10369;
    const Integer Q = 1152921504606538651ll;
    const Integer rootOfUnity = 199765477452278917ll;

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_ArbFFT12289(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[12289] = true;

    const usint N = 12289;
    const Integer Q = 1152921504606625421ll;
    const Integer rootOfUnity = 723240097479298971ll;

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_ArbFFT17497(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    primecyc::RaderFFTNat<Vector>::m_enabled[17497] = true;

    const usint N = 17497;
    const Integer Q = 1152921504606191953ll;
    const Integer rootOfUnity = 55945261116891098ll;

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N - 1);

    for (auto _ : state) {
        a = ChineseRemainderTransformArb<Vector>().ForwardTransform(a, rootOfUnity, 0, 0, N);
    }
}

static void BM_FFT23Vec(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 768;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);

    for (auto _ : state) {
        Vector b(N, Q);
        primecyc::RaderFFTNat<Vector>().ForwardFFTBase2n3(a, rootOfUnity, &b);
    }
}

static void BM_FFT23StdVec(benchmark::State &state) {
    using Vector = lbcrypto::NativeVector;
    using Integer = typename Vector::Integer;

    const usint N = 768;
    const Integer Q = lbcrypto::LastPrime<Integer>(59, N);
    const Integer rootOfUnity = lbcrypto::RootOfUnity(N, Q);

    Vector a = lbcrypto::DiscreteUniformGeneratorImpl<Vector>(Q).GenerateVector(N);
    std::vector<Integer> aStdVec(N);
    for (usint i = 0; i < N; i++) {
        aStdVec[i] = a[i];
    }

    for (auto _ : state) {
        std::vector<Integer> bStdVec(N);
        primecyc::RaderFFTNat<Vector>().ForwardFFTBase2n3StdVec(aStdVec, Q, rootOfUnity, bStdVec);
    }
}

BENCHMARK(BM_PowerOf2FFT512);
BENCHMARK(BM_PowerOf2FFT1024);
BENCHMARK(BM_PowerOf2FFT8192);
BENCHMARK(BM_PowerOf2FFT16384);

BENCHMARK(BM_ArbFFT769);
BENCHMARK(BM_ArbFFT1153);
BENCHMARK(BM_ArbFFT3889);
BENCHMARK(BM_ArbFFT10369);
BENCHMARK(BM_ArbFFT12289);
BENCHMARK(BM_ArbFFT17497);

BENCHMARK(BM_FFT23Vec);
BENCHMARK(BM_FFT23StdVec);

BENCHMARK_MAIN();