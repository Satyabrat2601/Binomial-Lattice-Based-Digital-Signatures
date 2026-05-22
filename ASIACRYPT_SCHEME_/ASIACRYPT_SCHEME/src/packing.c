
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <string.h>

/* ================================================================ Public Key */

void pack_pk(uint8_t pk[DUET_CRYPTO_PUBLICKEYBYTES],
             const uint8_t rho[DUET_SEEDBYTES], const polyveck *b) {
    memcpy(pk, rho, DUET_SEEDBYTES);
    for (unsigned i = 0; i < DUET_K; ++i)
        poly_pack_q(pk + DUET_SEEDBYTES + i*DUET_POLYQ_PACKEDBYTES, &b->vec[i]);
}

void unpack_pk(uint8_t rho[DUET_SEEDBYTES], polyveck *b,
               const uint8_t pk[DUET_CRYPTO_PUBLICKEYBYTES]) {
    memcpy(rho, pk, DUET_SEEDBYTES);
    for (unsigned i = 0; i < DUET_K; ++i)
        poly_unpack_q(&b->vec[i], pk + DUET_SEEDBYTES + i*DUET_POLYQ_PACKEDBYTES);
}

/* ================================================================ Secret Key */

void pack_sk(uint8_t sk[DUET_CRYPTO_SECRETKEYBYTES],
             const uint8_t rho[DUET_SEEDBYTES],
             const uint8_t key[DUET_SEEDBYTES],
             const uint8_t tr[DUET_SEEDBYTES],
             const polyvecell *s1, const polyveck *s2) {
    uint8_t *p = sk;
    memcpy(p, rho, DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    memcpy(p, key, DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    memcpy(p, tr,  DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    for (unsigned i = 0; i < DUET_ELL; ++i) { poly_pack_eta(p, &s1->vec[i]); p += DUET_POLYETA_PACKEDBYTES; }
    for (unsigned i = 0; i < DUET_K;   ++i) { poly_pack_eta(p, &s2->vec[i]); p += DUET_POLYETA_PACKEDBYTES; }
}

void unpack_sk(uint8_t rho[DUET_SEEDBYTES], uint8_t key[DUET_SEEDBYTES],
               uint8_t tr[DUET_SEEDBYTES], polyvecell *s1, polyveck *s2,
               const uint8_t sk[DUET_CRYPTO_SECRETKEYBYTES]) {
    const uint8_t *p = sk;
    memcpy(rho, p, DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    memcpy(key, p, DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    memcpy(tr,  p, DUET_SEEDBYTES); p += DUET_SEEDBYTES;
    for (unsigned i = 0; i < DUET_ELL; ++i) { poly_unpack_eta(&s1->vec[i], p); p += DUET_POLYETA_PACKEDBYTES; }
    for (unsigned i = 0; i < DUET_K;   ++i) { poly_unpack_eta(&s2->vec[i], p); p += DUET_POLYETA_PACKEDBYTES; }
}

/* ================================================================ Signature */

void pack_sig(uint8_t sig[DUET_CRYPTO_BYTES],
              const polyvecL *z1, const polyveck *h,
              const uint8_t c_seed[DUET_SEEDBYTES]) {
    uint8_t *p = sig;
    /* z1: (ell+1) polys at B_BITS=10 bits/coeff */
    for (unsigned i = 0; i <= DUET_ELL; ++i) { poly_pack_z(p, &z1->vec[i]); p += DUET_POLYZ_PACKEDBYTES; }
    /* h: k polys at 8 bits/coeff (HighBits correction values in [0,H_RANGE)) */
    for (unsigned i = 0; i < DUET_K; ++i) {
        for (unsigned r = 0; r < DUET_N; ++r)
            p[r] = (uint8_t)(h->vec[i].coeffs[r] & 0xFF);
        p += DUET_POLYH_PACKEDBYTES;
    }
    /* c_seed: 32 bytes */
    memcpy(p, c_seed, DUET_POLYC_PACKEDBYTES);
}

int unpack_sig(polyvecL *z1, polyveck *h, uint8_t c_seed[DUET_SEEDBYTES],
               const uint8_t sig[DUET_CRYPTO_BYTES]) {
    const uint8_t *p = sig;
    for (unsigned i = 0; i <= DUET_ELL; ++i) { poly_unpack_z(&z1->vec[i], p); p += DUET_POLYZ_PACKEDBYTES; }
    for (unsigned i = 0; i < DUET_K; ++i) {
        for (unsigned r = 0; r < DUET_N; ++r)
            h->vec[i].coeffs[r] = (int32_t)(p[r]);
        p += DUET_POLYH_PACKEDBYTES;
    }
    memcpy(c_seed, p, DUET_POLYC_PACKEDBYTES);
    return 0;
}
