
#include "polyvec.h"
#include "poly.h"
#include "params.h"
#include "reduce.h"
#include "symmetric.h"
#include <stdint.h>
#include <string.h>

/* ================================================================
   polyveck
   ================================================================ */

void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}
void polyveck_freeze(polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_freeze(&v->vec[i]);
}
void polyveck_freeze2q(polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_freeze2q(&v->vec[i]);
}
void polyveck_reduce2q(polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_reduce2q(&v->vec[i]);
}
void polyveck_ntt(polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_ntt(&v->vec[i]);
}
void polyveck_invntt_tomont(polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_invntt_tomont(&v->vec[i]);
}
void polyveck_highbits(polyveck *h, const polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) poly_highbits(&h->vec[i], &v->vec[i]);
}
int polyveck_chknorm(const polyveck *v, int32_t bound) {
    for (unsigned i = 0; i < DUET_K; ++i)
        if (poly_chknorm(&v->vec[i], bound)) return 1;
    return 0;
}
void polyveck_expand_eta(polyveck *v, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    for (unsigned i = 0; i < DUET_K; ++i)
        poly_uniform_eta(&v->vec[i], seed, nonce + (uint16_t)i);
}
void polyveck_pack_hi(uint8_t *buf, const polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i) {
        poly hi;
        poly_highbits(&hi, &v->vec[i]);
        for (unsigned r = 0; r < DUET_N; ++r)
            buf[i * DUET_N + r] = (uint8_t)(hi.coeffs[r] & 0xFF);
    }
}

void polyveck_unpack_hi(polyveck *v, const uint8_t *buf) {
    for (unsigned i = 0; i < DUET_K; ++i)
        for (unsigned r = 0; r < DUET_N; ++r)
            v->vec[i].coeffs[r] = (int32_t)(buf[i * DUET_N + r]);
}

void polyveck_pw_acc_mont(polyveck *t, const polyveck *v, const poly *c) {
    for (unsigned i = 0; i < DUET_K; ++i)
        poly_pointwise_montgomery(&t->vec[i], &v->vec[i], c);
}
void polyveck_scalar_mul_q(polyveck *out, const polyveck *v) {
    for (unsigned i = 0; i < DUET_K; ++i)
        poly_scalar_mul_q(&out->vec[i], &v->vec[i], DUET_Q);
}

/*
 * polyveck_sub_2z2 -- w' = w - 2*z2 mod 2q  (Fig.2 step 10)
 * w in [0,2q). z2 coefficients are signed (from Split); lift first.
 */
void polyveck_sub_2z2(polyveck *wp, const polyveck *w, const polyveck *z2) {
    int32_t dq = (int32_t)DUET_DQ;
    for (unsigned i = 0; i < DUET_K; ++i) {
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t wi  = w->vec[i].coeffs[r];
            int32_t z2i = z2->vec[i].coeffs[r];
            if (z2i < 0) z2i += dq;
            int32_t val = wi - 2 * z2i;
            val = ((val % dq) + dq) % dq;
            wp->vec[i].coeffs[r] = val;
        }
    }
}

/*
 * polyveck_make_hint -- h = MakeHint(w, w')
 * h[i][r] = (HighBits(w[i][r]) - HighBits(w'[i][r])) mod H_RANGE
 * Stored as small integer in [0, H_RANGE).
 * Matches Python estimator MakeHint exactly.
 */
void polyveck_make_hint(polyveck *h, const polyveck *w, const polyveck *wp) {
    for (unsigned i = 0; i < DUET_K; ++i) {
        poly hi_w, hi_wp;
        poly_highbits(&hi_w,  &w->vec[i]);
        poly_highbits(&hi_wp, &wp->vec[i]);
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t diff = hi_w.coeffs[r] - hi_wp.coeffs[r];
            diff = ((diff % (int32_t)DUET_H_RANGE) + (int32_t)DUET_H_RANGE)
                   % (int32_t)DUET_H_RANGE;
            h->vec[i].coeffs[r] = diff;
        }
    }
}

/* ================================================================
   polyvecell
   ================================================================ */

void polyvecell_add(polyvecell *w, const polyvecell *u, const polyvecell *v) {
    for (unsigned i = 0; i < DUET_ELL; ++i) poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}
void polyvecell_ntt(polyvecell *v) {
    for (unsigned i = 0; i < DUET_ELL; ++i) poly_ntt(&v->vec[i]);
}
void polyvecell_invntt_tomont(polyvecell *v) {
    for (unsigned i = 0; i < DUET_ELL; ++i) poly_invntt_tomont(&v->vec[i]);
}
int polyvecell_chknorm(const polyvecell *v, int32_t bound) {
    for (unsigned i = 0; i < DUET_ELL; ++i)
        if (poly_chknorm(&v->vec[i], bound)) return 1;
    return 0;
}
void polyvecell_expand_eta(polyvecell *v, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce) {
    for (unsigned i = 0; i < DUET_ELL; ++i)
        poly_uniform_eta(&v->vec[i], seed, nonce + (uint16_t)i);
}
void polyvecell_pw_acc_mont(polyvecell *t, const polyvecell *v, const poly *c) {
    for (unsigned i = 0; i < DUET_ELL; ++i)
        poly_pointwise_montgomery(&t->vec[i], &v->vec[i], c);
}

/* ================================================================
   polyvecL
   ================================================================ */

void polyvecL_add(polyvecL *w, const polyvecL *u, const polyvecL *v) {
    for (unsigned i = 0; i < DUET_L; ++i) poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/* Infinity-norm check: z coefficients are signed, no mod-q centering. */
int polyvecL_chknorm(const polyvecL *v, int32_t bound) {
    for (unsigned i = 0; i < DUET_L; ++i)
        if (poly_chknorm(&v->vec[i], bound)) return 1;
    return 0;
}

/* CBD_B validity check: |coeff| <= bound for all L*N coefficients. */
int polyvecL_is_in_cbd_bound(const polyvecL *v, int32_t bound) {
    for (unsigned i = 0; i < DUET_L; ++i)
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t c = v->vec[i].coeffs[r];
            if (c < 0) c = -c;
            if (c > bound) return 0;
        }
    return 1;
}

/* Sample y <- CBD_{gamma1}^L from seed (one poly per nonce). */
void polyvecL_expand_cbd_gamma1(polyvecL *y, const uint8_t seed[DUET_SEEDBYTES]) {
    for (unsigned i = 0; i < DUET_L; ++i)
        poly_cbd_gamma1(&y->vec[i], seed, (uint16_t)i);
}

void polyvecL_mul_xiS(polyvecL       *z,
                      const polyvecL *y,
                      const poly     *xi,
                      const polyvecL *S)
{
    /* Lead poly: z[0] = y[0] + xi (ring multiply by j=1 is identity) */
    for (unsigned r = 0; r < DUET_N; ++r)
        z->vec[0].coeffs[r] = y->vec[0].coeffs[r] + xi->coeffs[r];

    /* Remaining polys: z[j] = y[j] + negacyclic_conv(xi, S[j]) */
    for (unsigned j = 1; j < DUET_L; ++j) {
        int64_t conv[DUET_N];
        memset(conv, 0, sizeof(conv));
        for (unsigned p = 0; p < DUET_N; ++p) {
            int32_t xi_p = xi->coeffs[p];
            if (xi_p == 0) continue;  /* xi is sparse */
            for (unsigned q2 = 0; q2 < DUET_N; ++q2) {
                int64_t prod = (int64_t)xi_p * S->vec[j].coeffs[q2];
                unsigned idx = p + q2;
                if (idx < DUET_N) conv[idx]         += prod;
                else              conv[idx - DUET_N] -= prod;   /* negacyclic wrap */
            }
        }
        for (unsigned r = 0; r < DUET_N; ++r)
            z->vec[j].coeffs[r] = y->vec[j].coeffs[r] + (int32_t)conv[r];
    }
}

/*
 * polyvecL_split -- Split(z) -> (z1, z2)
 * z1 = z[0..ell]   (ell+1 polys, transmitted in signature)
 * z2 = z[ell+1..L-1] (k polys, used for MakeHint only, not transmitted)
 */
void polyvecL_split(polyvecL *z1, polyveck *z2, const polyvecL *z) {
    for (unsigned i = 0; i <= DUET_ELL; ++i)
        z1->vec[i] = z->vec[i];
    /* Zero unused polyvecL slots */
    for (unsigned i = DUET_ELL + 1; i < DUET_L; ++i)
        memset(&z1->vec[i], 0, sizeof(poly));
    for (unsigned i = 0; i < DUET_K; ++i)
        z2->vec[i] = z->vec[DUET_ELL + 1 + i];
}

