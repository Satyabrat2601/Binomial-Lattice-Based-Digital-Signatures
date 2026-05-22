
#ifndef DUET_POLYVEC_H
#define DUET_POLYVEC_H

#include "params.h"
#include "poly.h"
#include <stdint.h>

/* ---- polyveck (length k) ---- */
typedef struct { poly vec[DUET_K]; } polyveck;

#define polyveck_add              DUET_NAMESPACE(polyveck_add)
#define polyveck_sub              DUET_NAMESPACE(polyveck_sub)
#define polyveck_freeze           DUET_NAMESPACE(polyveck_freeze)
#define polyveck_freeze2q         DUET_NAMESPACE(polyveck_freeze2q)
#define polyveck_reduce2q         DUET_NAMESPACE(polyveck_reduce2q)
#define polyveck_ntt              DUET_NAMESPACE(polyveck_ntt)
#define polyveck_invntt_tomont    DUET_NAMESPACE(polyveck_invntt_tomont)
#define polyveck_highbits         DUET_NAMESPACE(polyveck_highbits)
#define polyveck_chknorm          DUET_NAMESPACE(polyveck_chknorm)
#define polyveck_expand_eta       DUET_NAMESPACE(polyveck_expand_eta)
#define polyveck_pw_acc_mont      DUET_NAMESPACE(polyveck_pw_acc_mont)
#define polyveck_scalar_mul_q     DUET_NAMESPACE(polyveck_scalar_mul_q)
#define polyveck_sub_2z2          DUET_NAMESPACE(polyveck_sub_2z2)
#define polyveck_make_hint        DUET_NAMESPACE(polyveck_make_hint)
#define polyveck_pack_hi          DUET_NAMESPACE(polyveck_pack_hi)
#define polyveck_unpack_hi        DUET_NAMESPACE(polyveck_unpack_hi)

void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v);
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v);
void polyveck_freeze(polyveck *v);
void polyveck_freeze2q(polyveck *v);
void polyveck_reduce2q(polyveck *v);
void polyveck_ntt(polyveck *v);
void polyveck_invntt_tomont(polyveck *v);
void polyveck_highbits(polyveck *h, const polyveck *v);
int  polyveck_chknorm(const polyveck *v, int32_t bound);
void polyveck_expand_eta(polyveck *v, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
void polyveck_pw_acc_mont(polyveck *t, const polyveck *v, const poly *c);
void polyveck_scalar_mul_q(polyveck *out, const polyveck *v);

/* w' = w - 2*z2 mod 2q */
void polyveck_sub_2z2(polyveck *wp, const polyveck *w, const polyveck *z2);

/* MakeHint: h[i][r] = (HighBits(w) - HighBits(w')) mod H_RANGE */
void polyveck_make_hint(polyveck *h, const polyveck *w, const polyveck *wp);

void polyveck_pack_hi(uint8_t *buf, const polyveck *v);
void polyveck_unpack_hi(polyveck *v, const uint8_t *buf);

/* ---- polyvecell (length ell) ---- */
typedef struct { poly vec[DUET_ELL]; } polyvecell;

#define polyvecell_add            DUET_NAMESPACE(polyvecell_add)
#define polyvecell_ntt            DUET_NAMESPACE(polyvecell_ntt)
#define polyvecell_invntt_tomont  DUET_NAMESPACE(polyvecell_invntt_tomont)
#define polyvecell_chknorm        DUET_NAMESPACE(polyvecell_chknorm)
#define polyvecell_expand_eta     DUET_NAMESPACE(polyvecell_expand_eta)
#define polyvecell_pw_acc_mont    DUET_NAMESPACE(polyvecell_pw_acc_mont)

void polyvecell_add(polyvecell *w, const polyvecell *u, const polyvecell *v);
void polyvecell_ntt(polyvecell *v);
void polyvecell_invntt_tomont(polyvecell *v);
int  polyvecell_chknorm(const polyvecell *v, int32_t bound);
void polyvecell_expand_eta(polyvecell *v, const uint8_t seed[DUET_SEEDBYTES], uint16_t nonce);
void polyvecell_pw_acc_mont(polyvecell *t, const polyvecell *v, const poly *c);

/* ---- polyvecL (length L = k+ell+1) ---- */
typedef struct { poly vec[DUET_L]; } polyvecL;

#define polyvecL_add                DUET_NAMESPACE(polyvecL_add)
#define polyvecL_chknorm            DUET_NAMESPACE(polyvecL_chknorm)
#define polyvecL_is_in_cbd_bound    DUET_NAMESPACE(polyvecL_is_in_cbd_bound)
#define polyvecL_expand_cbd_gamma1  DUET_NAMESPACE(polyvecL_expand_cbd_gamma1)
#define polyvecL_mul_xiS            DUET_NAMESPACE(polyvecL_mul_xiS)
#define polyvecL_split              DUET_NAMESPACE(polyvecL_split)

void polyvecL_add(polyvecL *w, const polyvecL *u, const polyvecL *v);
int  polyvecL_chknorm(const polyvecL *v, int32_t bound);
int  polyvecL_is_in_cbd_bound(const polyvecL *v, int32_t bound);
void polyvecL_expand_cbd_gamma1(polyvecL *y, const uint8_t seed[DUET_SEEDBYTES]);
void polyvecL_mul_xiS(polyvecL *z, const polyvecL *y, const poly *xi, const polyvecL *S);
void polyvecL_split(polyvecL *z1, polyveck *z2, const polyvecL *z);

/* xi sampling from challenge */
#define poly_sample_xi  DUET_NAMESPACE(poly_sample_xi)
void poly_sample_xi(poly *xi, const poly *c, const uint8_t seed[DUET_SEEDBYTES]);

#endif /* DUET_POLYVEC_H */
