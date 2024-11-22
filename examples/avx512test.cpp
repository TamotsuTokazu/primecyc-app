#include <iostream>
#include <immintrin.h>  // For intrinsic functions
#include <cpuid.h>      // To use __get_cpuid

bool is_avx512_supported() {
    unsigned int eax, ebx, ecx, edx;

    // Query for extended CPU features (Leaf 7, Sub-leaf 0)
    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        // Check for the AVX-512 Foundation bit (bit 16 of EBX)
        if (ebx & (1 << 16)) {
            return true;  // AVX-512 is supported
        }
    }

    return false;  // AVX-512 is not supported
}

int main() {
    if (is_avx512_supported()) {
        std::cout << "AVX-512 is supported on this CPU." << std::endl;
    } else {
        std::cout << "AVX-512 is NOT supported on this CPU." << std::endl;
    }
    return 0;
}