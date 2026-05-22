
#ifndef DUET_POLY_H
#define DUET_POLY_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef struct { int32_t coeffs[DUET_N]; } poly;

/* ---- Arithmetic ---- */
#define poly_add                  DUET_NAMESPACE(poly_add)
#define poly_sub                  DUET_NAMESPACE(poly_sub)
#define poly_pointwise_montgomery DUET_NAMESPACE(poly_pointwise_montgomery)
#define poly_scalar_mul_q         DUET_NAMESPACE(poly_scalar_mul_q)
#define poly_tomont               DUET_NAMESPACE(poly_tomont)
#define poly_frommont             DUET_NAMESPACE(poly_frommont)

void poly_add(poly *c, const poly *a, const poly *b);
void poly_sub(poly *c, const poly *a, const poly *b);
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);
void poly_scalar_mul_q(poly *c, const poly *a, int32_t s);
void poly_tomont(poly *a);
void poly_frommont(poly *a);

/* ---- Reduction ---- */
/* poly_reduce: Barrett reduce each coeff to (-q, q) for intermediate NTT results */
#define poly_reduce     DUET_NAMESPACE(poly_reduce)
/* poly_reduce2q: reduce each coeff to [0, 2q) for the 2q-domain (A_hat operations) */
#define poly_reduce2q   DUET_NAMESPACE(poly_reduce2q)
#define poly_freeze    DUET_NAMESPACE(poly_freeze)
#define poly_freeze2q  DUET_NAMESPACE(poly_freeze2q)
#define poly_caddq     DUET_NAMESPACE(poly_caddq)

void poly_reduce(poly *a);
void poly_reduce2q(poly *a);
void poly_freeze(poly *a);
void poly_freeze2q(poly *a);
void poly_caddq(poly *a);

/* ---- NTT ---- */
#define poly_ntt           DUET_NAMESPACE(poly_ntt)
#define poly_invntt_tomont DUET_NAMESPACE(poly_invntt_tomont)

void poly_ntt(poly *a);
void poly_invntt_tomont(poly *a);

/* ---- Decomposition ---- */
#define poly_highbits  DUET_NAMESPACE(poly_highbits)
#define poly_lowbits   DUET_NAMESPACE(poly_lowbits)
#define poly_usehint   DUET_NAMESPACE(poly_usehint)

void poly_highbits(poly *hi, const poly *a);
void poly_lowbits(poly *lo, const poly *a);
void poly_usehint(poly *out, const poly *r, const poly *h);

/* ---- Sampling ---- */
#define poly_uniform        DUET_NAMESPACE(poly_uniform)
#define poly_uniform_eta    DUET_NAMESPACE(poly_uniform_eta)
#define poly_uniform_gamma1 DUET_NAMESPACE(poly_uniform_gamma1)
#define poly_cbd_gamma1     DUET_NAMESPACE(poly_cbd_gamma1)
#define poly_challenge      DUET_NAMESPACE(poly_challenge)
#define poly_sample_xi      DUET_NAMESPACE(poly_sample_xi)

void poly_uniform(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
void poly_uniform_eta(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
void poly_uniform_gamma1(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
void poly_cbd_gamma1(poly *a, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
/* Weight-tau {0,+1} challenge from SHAKE-256(w1 || w0 || mu) */
void poly_challenge(poly *c,
                    const uint8_t *w_packed, size_t w_packed_len,
                    const uint8_t mu_bar[DUET_SEEDBYTES]);
/* Xi sampler: xi_j=0 if c_j=0, else +-1 via u~CBD_1, d~Unif{0,1} */
void poly_sample_xi(poly *xi, const poly *c, const uint8_t seed[DUET_SEEDBYTES]);

/* ---- Norm ---- */
#define poly_chknorm  DUET_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t bound);

/* ---- Packing ---- */
#define poly_pack_q    DUET_NAMESPACE(poly_pack_q)
#define poly_unpack_q  DUET_NAMESPACE(poly_unpack_q)
#define poly_pack_eta  DUET_NAMESPACE(poly_pack_eta)
#define poly_unpack_eta DUET_NAMESPACE(poly_unpack_eta)
#define poly_pack_z    DUET_NAMESPACE(poly_pack_z)
#define poly_unpack_z  DUET_NAMESPACE(poly_unpack_z)

void poly_pack_q(uint8_t buf[DUET_POLYQ_PACKEDBYTES], const poly *a);
void poly_unpack_q(poly *a, const uint8_t buf[DUET_POLYQ_PACKEDBYTES]);
void poly_pack_eta(uint8_t buf[DUET_POLYETA_PACKEDBYTES], const poly *a);
void poly_unpack_eta(poly *a, const uint8_t buf[DUET_POLYETA_PACKEDBYTES]);
/* z1 coeffs in [-B,B] packed at B_BITS=10 bits, centered on B */
void poly_pack_z(uint8_t buf[DUET_POLYZ_PACKEDBYTES], const poly *a);
void poly_unpack_z(poly *a, const uint8_t buf[DUET_POLYZ_PACKEDBYTES]);

#endif /* !DUET_POLY_H */
