
#include "fips202.h"
#include "params.h"
#include "symmetric.h"
#include <stdint.h>

void shake128_stream_init(keccak_state *st,
                          const uint8_t seed[DUET_SEEDBYTES],
                          uint16_t nonce) {
    uint8_t t[2] = {(uint8_t)nonce, (uint8_t)(nonce >> 8)};
    shake128_init(st);
    shake128_absorb(st, seed, DUET_SEEDBYTES);
    shake128_absorb(st, t, 2);
    shake128_finalize(st);
}

void shake256_stream_init(keccak_state *st,
                          const uint8_t seed[DUET_SEEDBYTES],
                          uint16_t nonce) {
    uint8_t t[2] = {(uint8_t)nonce, (uint8_t)(nonce >> 8)};
    shake256_init(st);
    shake256_absorb(st, seed, DUET_SEEDBYTES);
    shake256_absorb(st, t, 2);
    shake256_finalize(st);
}

void shake256_absorb_twice(keccak_state *st,
                            const uint8_t *in1, size_t in1len,
                            const uint8_t *in2, size_t in2len) {
    shake256_init(st);
    shake256_absorb(st, in1, in1len);
    shake256_absorb(st, in2, in2len);
    shake256_finalize(st);
}

void shake256_absorb_thrice(keccak_state *st,
                             const uint8_t *in1, size_t in1len,
                             const uint8_t *in2, size_t in2len,
                             const uint8_t *in3, size_t in3len) {
    shake256_init(st);
    shake256_absorb(st, in1, in1len);
    shake256_absorb(st, in2, in2len);
    shake256_absorb(st, in3, in3len);
    shake256_finalize(st);
}
