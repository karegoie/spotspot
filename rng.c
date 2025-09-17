#include "rng.h"
#include <stdint.h>

/* MT19937 constants */
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908B0DFUL
#define MT_UPPER_MASK 0x80000000UL
#define MT_LOWER_MASK 0x7FFFFFFFUL

/* State */
static uint32_t mt[MT_N];
static int mt_index = MT_N + 1; /* mt_index == N+1 means state not initialized */

static void mt_seed(uint32_t s) {
    mt[0] = s;
    for (mt_index = 1; mt_index < MT_N; mt_index++) {
        mt[mt_index] = (1812433253UL * (mt[mt_index - 1] ^ (mt[mt_index - 1] >> 30)) + (uint32_t)mt_index);
    }
}

static void mt_twist(void) {
    static const uint32_t mag01[2] = {0x0UL, MT_MATRIX_A};
    uint32_t y;
    int kk;

    for (kk = 0; kk < MT_N - MT_M; kk++) {
        y = (mt[kk] & MT_UPPER_MASK) | (mt[kk + 1] & MT_LOWER_MASK);
        mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < MT_N - 1; kk++) {
        y = (mt[kk] & MT_UPPER_MASK) | (mt[kk + 1] & MT_LOWER_MASK);
        mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[MT_N - 1] & MT_UPPER_MASK) | (mt[0] & MT_LOWER_MASK);
    mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mt_index = 0;
}

/* Extract a tempered 32-bit unsigned int */
static uint32_t mt_rand_uint32(void) {
    uint32_t y;

    if (mt_index >= MT_N) {
        if (mt_index == MT_N + 1) {
            /* default seed */
            mt_seed(5489UL);
        } else {
            mt_twist();
        }
    }

    y = mt[mt_index++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9D2C5680UL;
    y ^= (y << 15) & 0xEFC60000UL;
    y ^= (y >> 18);

    return y;
}

void init_rng(unsigned long seed) {
    mt_seed((uint32_t)seed);
}

/* Return double in [0,1) with 53-bit resolution using two 32-bit draws */
double random_double(void) {
    /* Construct a 53-bit mantissa: (a>>5)*2^26 + (b>>6) over 2^53 */
    uint32_t a = mt_rand_uint32();
    uint32_t b = mt_rand_uint32();
    uint64_t x = ((uint64_t)(a >> 5) << 26) | (uint64_t)(b >> 6);
    return (double)x / 9007199254740992.0; /* 2^53 */
}
