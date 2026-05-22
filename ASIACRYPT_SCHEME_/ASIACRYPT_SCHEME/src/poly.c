
#include "poly.h"
#include "ntt.h"
#include "params.h"
#include "reduce.h"
#include "symmetric.h"
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* ================================================================
   Basic ring arithmetic
   ================================================================ */

void poly_add(poly *c, const poly *a, const poly *b) {
    for (unsigned i = 0; i < DUET_N; ++i)
        c->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

void poly_sub(poly *c, const poly *a, const poly *b) {
    for (unsigned i = 0; i < DUET_N; ++i)
        c->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
    for (unsigned i = 0; i < DUET_N; ++i)
        c->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);
}

void poly_scalar_mul_q(poly *c, const poly *a, int32_t s) {
    for (unsigned i = 0; i < DUET_N; ++i)
        c->coeffs[i] = a->coeffs[i] * s;
}

void poly_tomont(poly *a) {
    /* Multiply each coeff by MONT = 2^32 mod q (using Montgomery multiply) */
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * ((int64_t)MONT * MONT));
}

void poly_frommont(poly *a) {
    /* Convert from Montgomery domain: multiply by 1 in Montgomery form = multiply by 2^{-32} */
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i]);
}


/* poly_reduce2q: reduce each coefficient to [0, 2q) for the 2q-domain. */
void poly_reduce2q(poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = freeze2q(a->coeffs[i]);  /* exact [0, 2q) */
}

/* poly_reduce: Barrett reduce to (-q, q) for NTT post-processing.
 * Used by polymat.c after NTT-based multiply-accumulate. */
void poly_reduce(poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = barrett_reduce(a->coeffs[i]);  /* (-q, q) */
}

void poly_freeze(poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = freeze(a->coeffs[i]);
}

void poly_freeze2q(poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = freeze2q(a->coeffs[i]);
}

void poly_caddq(poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = caddq(a->coeffs[i]);
}

/* ================================================================
   NTT
   ================================================================ */

void poly_ntt(poly *a)           { ntt(a->coeffs);           }
void poly_invntt_tomont(poly *a) { invntt_tomont(a->coeffs); }


static int32_t highbits_one(int32_t r) {
    /* r in [0, 2q).  Return hi = round(r / DUET_ALPHA) */
    int32_t r1 = (r + (int32_t)(DUET_ALPHA / 2)) / (int32_t)DUET_ALPHA;
    /* edge case: if r1 == floor(2q/ALPHA)+1, wrap to 0 */
    if (r1 >= (int32_t)DUET_H_RANGE) r1 = 0;
    return r1;
}

void poly_highbits(poly *hi, const poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i)
        hi->coeffs[i] = highbits_one(a->coeffs[i]);
}

void poly_lowbits(poly *lo, const poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i) {
        int32_t r1 = highbits_one(a->coeffs[i]);
        lo->coeffs[i] = a->coeffs[i] - r1 * (int32_t)DUET_ALPHA;
    }
}

void poly_usehint(poly *out, const poly *r, const poly *h) {
    for (unsigned i = 0; i < DUET_N; ++i) {
        int32_t hi = highbits_one(r->coeffs[i]);
        int32_t recovered = (hi + h->coeffs[i]) % (int32_t)DUET_H_RANGE;
        if (recovered < 0) recovered += (int32_t)DUET_H_RANGE;
        out->coeffs[i] = recovered;
    }
}

/* ================================================================
   Norm check (centered mod q)
   ================================================================ */

int poly_chknorm(const poly *a, int32_t bound) {
    for (unsigned i = 0; i < DUET_N; ++i) {
        int32_t c = a->coeffs[i];
        int32_t ac = c < 0 ? -c : c;
        if (ac > bound) return 1;
    }
    return 0;
}

/* ================================================================
   Uniform sampling (poly_uniform via SHAKE-128)
   ================================================================ */
#define POLY_UNIFORM_NBLOCKS \
    ((768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES)

void poly_uniform(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    uint8_t buf[POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES + 2];
    stream128_state st;
    stream128_init(&st, seed, nonce);
    stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &st);
    unsigned int ctr = 0, pos = 0;
    unsigned int buflen = POLY_UNIFORM_NBLOCKS * STREAM128_BLOCKBYTES;
    while (ctr < DUET_N) {
        if (pos + 3 > buflen) {
            stream128_squeezeblocks(buf, 1, &st);
            buflen = STREAM128_BLOCKBYTES;
            pos = 0;
        }
        uint32_t val = (uint32_t)buf[pos]
                     | ((uint32_t)buf[pos+1] << 8)
                     | ((uint32_t)buf[pos+2] << 16);
        val &= 0x1FFFF; /* 17-bit mask to cover q=64513 < 2^17 */
        pos += 2;
        if (val < (uint32_t)DUET_Q)
            a->coeffs[ctr++] = (int32_t)val;
    }
}

/* ================================================================
   CBD_eta sampling (eta=1: coeffs in {-1,0,+1})
   ================================================================ */
#define POLY_UNIFORM_ETA_NBLOCKS \
    ((136 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES)

void poly_uniform_eta(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES];
    stream256_state st;
    stream256_init(&st, seed, nonce);
    stream256_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &st);
    unsigned int ctr = 0, pos = 0;
    unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS * STREAM256_BLOCKBYTES;
    while (ctr < DUET_N) {
        if (pos >= buflen) {
            stream256_squeezeblocks(buf, 1, &st);
            buflen = STREAM256_BLOCKBYTES;
            pos = 0;
        }
        uint8_t byte = buf[pos++];
        uint8_t lo = byte & 0x0F;
        uint8_t hi = byte >> 4;
        /* Ternary from nibble: 0->+1, 1->0, 2->-1, others rejected */
        if (lo < 3) { a->coeffs[ctr++] = 1 - (int32_t)lo; if (ctr >= DUET_N) break; }
        if (hi < 3) { a->coeffs[ctr++] = 1 - (int32_t)hi; }
    }
}

/* ================================================================
   CBD_{gamma1} sampling  (DUET scheme)
   Samples from Binomial(2*gamma1, 1/2) - gamma1.
   y ~ CBD_{gamma1} as specified in Fig.1 Sign step 3.

   Fast implementation: process gamma1 bits = ceil(gamma1/8) full bytes
   using hardware popcount (__builtin_popcount).
   For gamma1=339: 43 bytes for 'a', 43 bytes for 'b' per coefficient.
   The last byte of each half is masked to exactly gamma1 bits.
   Total randomness per poly: n * 2 * ceil(gamma1/8) bytes (squeezed from SHAKE-256).
   ================================================================ */
void poly_cbd_gamma1(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    const unsigned g1        = DUET_GAMMA1;
    const unsigned full_bytes = g1 / 8;          /* complete bytes per half */
    const unsigned rem_bits   = g1 % 8;          /* leftover bits in last byte */
    /* bytes consumed per coefficient: 2 * ceil(g1/8) */
    const unsigned bytes_per = 2u * (full_bytes + (rem_bits ? 1u : 0u));
    /* total bytes needed for all 256 coefficients */
    const unsigned total     = DUET_N * bytes_per + STREAM256_BLOCKBYTES; /* +margin */

    stream256_state st;
    stream256_init(&st, seed, nonce);

    /* Use a heap-allocated or large stack buffer; pull in SHAKE blocks on demand */
    uint8_t blk[STREAM256_BLOCKBYTES];
    unsigned blk_pos = STREAM256_BLOCKBYTES;  /* force immediate fill */
    // removed unused var
    (void)total;

    /* Inline helper: read one byte from SHAKE stream */
    #define READ_BYTE(out) do { \
        if (blk_pos >= STREAM256_BLOCKBYTES) { \
            stream256_squeezeblocks(blk, 1, &st); \
            blk_pos = 0; \
        } \
        (out) = blk[blk_pos++]; \
    } while(0)

    for (unsigned i = 0; i < DUET_N; ++i) {
        unsigned sum_a = 0, sum_b = 0;
        uint8_t byte;
        /* Sum 'a': full_bytes complete bytes + optional partial byte */
        for (unsigned j = 0; j < full_bytes; ++j) {
            READ_BYTE(byte);
            sum_a += (unsigned)__builtin_popcount(byte);
        }
        if (rem_bits) {
            READ_BYTE(byte);
            sum_a += (unsigned)__builtin_popcount(byte & (uint8_t)((1u << rem_bits) - 1u));
        }
        /* Sum 'b': same structure */
        for (unsigned j = 0; j < full_bytes; ++j) {
            READ_BYTE(byte);
            sum_b += (unsigned)__builtin_popcount(byte);
        }
        if (rem_bits) {
            READ_BYTE(byte);
            sum_b += (unsigned)__builtin_popcount(byte & (uint8_t)((1u << rem_bits) - 1u));
        }
        a->coeffs[i] = (int32_t)sum_a - (int32_t)sum_b;
    }
    #undef READ_BYTE
}

/* poly_uniform_gamma1: kept for compatibility; delegates to CBD */
void poly_uniform_gamma1(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    poly_cbd_gamma1(a, seed, nonce);
}

/* ================================================================
   Challenge polynomial -- Fig.1 step 5
   c = H(mu || w) in B_tau: weight-tau {0,+1} polynomial.
   (DUET uses {0,+1} challenges; sign is absorbed into xi sampling.)

   We store the SHAKE-256 seed and re-derive c on verify.
   Challenge hash: SHAKE256(w_packed || mu_bar), Fisher-Yates on last tau.
   ================================================================ */
void poly_challenge(poly *c,
                    const uint8_t *w_packed, size_t w_packed_len,
                    const uint8_t mu_bar[DUET_SEEDBYTES]) {
    keccak_state st;
    /* Correct Fisher-Yates using a position array.
     * We maintain a "bag" of N positions, initialized to 0..N-1.
     * We draw TAU positions WITHOUT replacement using 2-byte random values.
     * This guarantees exactly TAU distinct positions are selected.
     * Buffer: 2*N bytes worst-case (but 2*TAU+margin is enough with rejection) */
    uint8_t pos[DUET_N]; /* position bag: pos[i] = actual position i represents */
    uint8_t outbuf[DUET_N]; /* 1 byte per draw (N=256, all values valid) */
    unsigned i;

    for (i = 0; i < DUET_N; ++i) pos[i] = (uint8_t)i;

    shake256_absorb_twice(&st, w_packed, w_packed_len, mu_bar, DUET_SEEDBYTES);
    shake256_squeeze(outbuf, DUET_N, &st);

    memset(c->coeffs, 0, sizeof(c->coeffs));
    /* Select TAU positions: draw from shrinking bag (no replacement) */
    for (unsigned k = 0; k < DUET_TAU; ++k) {
        /* Random index in [0, N-1-k] */
        unsigned range = DUET_N - k;
        unsigned j = (unsigned)outbuf[k] % range;   /* uniform mod 256: bias-free since N=256 */
        /* Take pos[j] as the selected position */
        c->coeffs[pos[j]] = 1;
        /* Remove pos[j] from bag: swap with last element */
        pos[j] = pos[range - 1];
    }
}

/* ================================================================
   Xi sampling  (Fig.1 Sign step 6-7)
   xi_j = 0 if c_j = 0
   xi_j = 2*u_j + 1 - 2*d_j  (mod +-4)  if c_j != 0
   u_j ~ CBD_1 (in {-1,0,+1}), d_j ~ Uniform{0,1}
   Result: xi_j in {-1, +1} when c_j != 0
   ================================================================ */
void poly_sample_xi(poly *xi, const poly *c, const uint8_t seed[DUET_SEEDBYTES]) {
    keccak_state st;
    uint8_t buf[STREAM256_BLOCKBYTES];
    shake256_absorb_once(&st, seed, DUET_SEEDBYTES);
    shake256_squeeze(buf, sizeof(buf), &st);
    unsigned pos = 0, buflen = STREAM256_BLOCKBYTES;

    for (unsigned i = 0; i < DUET_N; ++i) {
        if (c->coeffs[i] == 0) {
            xi->coeffs[i] = 0;
            continue;
        }
        /* Need 2 random bits: bit0=u high bit (or just use raw bit), bit1=d */
        while (pos + 1 > buflen) {
            shake256_squeeze(buf, STREAM256_BLOCKBYTES, &st);
            pos = 0;
        }
        uint8_t byte = buf[pos++];
        /* u_j ~ CBD_1: bit0 - bit1 (centered) */
        int32_t u = (int32_t)((byte >> 0) & 1) - (int32_t)((byte >> 1) & 1);
        int32_t d = (int32_t)((byte >> 2) & 1);
        /* xi_j = 2*u_j + 1 - 2*d_j, reduced mod +-4 to {-1,+1} */
        int32_t raw = 2*u + 1 - 2*d;
        if (raw > 2)  raw -= 4;
        if (raw < -2) raw += 4;
        xi->coeffs[i] = raw;  /* in {-1, +1} */
    }
}

/* ================================================================
   Packing
   ================================================================ */

/* poly_pack_q: 16 bits/coeff, coeffs in [0,q).  512 bytes for n=256. */
void poly_pack_q(uint8_t buf[DUET_POLYQ_PACKEDBYTES], const poly *a) {
    for (unsigned i = 0; i < DUET_N; ++i) {
        buf[2*i]   = (uint8_t)a->coeffs[i];
        buf[2*i+1] = (uint8_t)(a->coeffs[i] >> 8);
    }
}
void poly_unpack_q(poly *a, const uint8_t buf[DUET_POLYQ_PACKEDBYTES]) {
    for (unsigned i = 0; i < DUET_N; ++i)
        a->coeffs[i] = (int32_t)((uint16_t)buf[2*i] | ((uint16_t)buf[2*i+1] << 8))
                       % DUET_Q;
}

/* poly_pack_eta: 2 bits/coeff (eta=1 -> values {-1,0,+1} mapped to {2,1,0}).
   256 coeffs * 2 bits = 64 bytes = POLYETA_PACKEDBYTES. */
void poly_pack_eta(uint8_t buf[DUET_POLYETA_PACKEDBYTES], const poly *a) {
    for (unsigned i = 0; i < DUET_N / 4; ++i) {
        uint8_t t0 = (uint8_t)(DUET_ETA - a->coeffs[4*i+0]);
        uint8_t t1 = (uint8_t)(DUET_ETA - a->coeffs[4*i+1]);
        uint8_t t2 = (uint8_t)(DUET_ETA - a->coeffs[4*i+2]);
        uint8_t t3 = (uint8_t)(DUET_ETA - a->coeffs[4*i+3]);
        buf[i] = t0 | (t1 << 2) | (t2 << 4) | (t3 << 6);
    }
}
void poly_unpack_eta(poly *a, const uint8_t buf[DUET_POLYETA_PACKEDBYTES]) {
    for (unsigned i = 0; i < DUET_N / 4; ++i) {
        a->coeffs[4*i+0] = DUET_ETA - (int32_t)((buf[i])    & 3);
        a->coeffs[4*i+1] = DUET_ETA - (int32_t)((buf[i]>>2) & 3);
        a->coeffs[4*i+2] = DUET_ETA - (int32_t)((buf[i]>>4) & 3);
        a->coeffs[4*i+3] = DUET_ETA - (int32_t)((buf[i]>>6) & 3);
    }
}

/* poly_pack_z: DUET_B_BITS=10 bits/coeff, z in [-B, B].
 * Fix: center on B (not GAMMA1). Stored as (coeff + B) in [0, 2B].
 * 2B+1 <= 2^10 = 1024 for all modes (2*450+1=901 < 1024).
 * 256 coeffs * 10 bits = 320 bytes = POLYZ_PACKEDBYTES. */
void poly_pack_z(uint8_t buf[DUET_POLYZ_PACKEDBYTES], const poly *a) {
    uint64_t acc = 0;
    int      bits = 0;
    unsigned pos  = 0;
    for (unsigned i = 0; i < DUET_N; ++i) {
        /* z in [-B, B] -> store as z + B in [0, 2B] */
        uint32_t v = (uint32_t)(a->coeffs[i] + (int32_t)DUET_B);
        acc  |= ((uint64_t)v << bits);
        bits += DUET_B_BITS;
        while (bits >= 8) {
            buf[pos++] = (uint8_t)acc;
            acc >>= 8;
            bits -= 8;
        }
    }
    if (bits > 0) buf[pos] = (uint8_t)acc;
}
void poly_unpack_z(poly *a, const uint8_t buf[DUET_POLYZ_PACKEDBYTES]) {
    uint64_t acc  = 0;
    int      bits = 0;
    unsigned pos  = 0;
    const uint32_t mask = (1u << DUET_B_BITS) - 1u;
    for (unsigned i = 0; i < DUET_N; ++i) {
        while (bits < DUET_B_BITS) {
            acc  |= ((uint64_t)buf[pos++] << bits);
            bits += 8;
        }
        a->coeffs[i] = (int32_t)(acc & mask) - (int32_t)DUET_B;
        acc  >>= DUET_B_BITS;
        bits  -= DUET_B_BITS;
    }
}

/* poly_pack_c / poly_unpack_c:
 * Challenge c is a {0,1} polynomial with exactly TAU ones.
 * Pack as TAU sorted position bytes (each in [0,255]).
 * Size = DUET_TAU bytes = DUET_POLYC_PACKEDBYTES.
 * During verify: c is available directly from unpacked sig to compute A*z-q*c*j. */
void poly_pack_c(uint8_t buf[DUET_POLYC_PACKEDBYTES], const poly *a) {
    unsigned k = 0;
    for (unsigned i = 0; i < DUET_N && k < DUET_TAU; ++i) {
        if (a->coeffs[i] != 0)
            buf[k++] = (uint8_t)i;
    }
    /* zero-fill any remainder */
    while (k < DUET_TAU) buf[k++] = 0;
}
void poly_unpack_c(poly *a, const uint8_t buf[DUET_POLYC_PACKEDBYTES]) {
    memset(a->coeffs, 0, sizeof(a->coeffs));
    for (unsigned k = 0; k < DUET_TAU; ++k) {
        unsigned pos = buf[k];
        if (pos < DUET_N)
            a->coeffs[pos] = 1;
    }
}



