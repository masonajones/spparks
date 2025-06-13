#ifndef XOSHIRO256_AVX2_H
#define XOSHIRO256_AVX2_H

#include <immintrin.h>
#include <stdint.h>

// SplitMix64 for seeding
inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

// Left rotation
inline __m256i rotl(__m256i x, int k) {
    return _mm256_or_si256(_mm256_slli_epi64(x, k), _mm256_srli_epi64(x, 64 - k));
}

// Xoshiro256++ AVX2
class Xoshiro256AVX2 {
private:
    __m256i s0, s1, s2, s3;
    
public:
    // Initialize with a scalar seed
    Xoshiro256AVX2(uint64_t seed) {
        uint64_t tmp[16];
        for (int i = 0; i < 16; i++)
            tmp[i] = splitmix64(seed);
        s0 = _mm256_set_epi64x(tmp[3], tmp[2], tmp[1], tmp[0]);
        s1 = _mm256_set_epi64x(tmp[7], tmp[6], tmp[5], tmp[4]);
        s2 = _mm256_set_epi64x(tmp[11], tmp[10], tmp[9], tmp[8]);
        s3 = _mm256_set_epi64x(tmp[15], tmp[14], tmp[13], tmp[12]);
    }

    // Generate 4 double precision values in [0, 1)

    __m256d next_batch() {
        __m256i result = _mm256_add_epi64(s0, s3);
    
        __m256i t = _mm256_slli_epi64(s1, 17);
        s2 = _mm256_xor_si256(s2, s0);
        s3 = _mm256_xor_si256(s3, s1);
        s1 = _mm256_xor_si256(s1, s2);
        s0 = _mm256_xor_si256(s0, s3);
        s2 = _mm256_xor_si256(s2, t);
        s3 = rotl(s3, 45);
    
        // Convert 64-bit int to double in [0,1)
        alignas(32) uint64_t raw[4];
        _mm256_store_si256((__m256i*)raw, result);
        return _mm256_set_pd(
            (raw[3] >> 11) * (1.0 / (1ULL << 53)),
            (raw[2] >> 11) * (1.0 / (1ULL << 53)),
            (raw[1] >> 11) * (1.0 / (1ULL << 53)),
            (raw[0] >> 11) * (1.0 / (1ULL << 53))
        );
    }
    
    void reset(uint64_t seed) {
        uint64_t tmp[16];
        for (int i = 0; i < 16; i++)
            tmp[i] = splitmix64(seed);
        s0 = _mm256_set_epi64x(tmp[3], tmp[2], tmp[1], tmp[0]);
        s1 = _mm256_set_epi64x(tmp[7], tmp[6], tmp[5], tmp[4]);
        s2 = _mm256_set_epi64x(tmp[11], tmp[10], tmp[9], tmp[8]);
        s3 = _mm256_set_epi64x(tmp[15], tmp[14], tmp[13], tmp[12]);
    }

};

#endif
