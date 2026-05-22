
#ifndef DUET_POLYMAT_H
#define DUET_POLYMAT_H

#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <stdint.h>

/* ---- Expand base matrix A from seed ---- */
#define polymat_expand_A  DUET_NAMESPACE(polymat_expand_A)
void polymat_expand_A(poly A[DUET_K][DUET_ELL], const uint8_t rho[DUET_SEEDBYTES]);

/* ---- Build A_hat from A and b ---- */
#define polymat_build_A_hat  DUET_NAMESPACE(polymat_build_A_hat)
void polymat_build_A_hat(polyvecL A_hat[DUET_K],
                         const poly A[DUET_K][DUET_ELL],
                         const polyveck *b);

/* ---- w = A_hat * z  (mod 2q) - ---- */
#define polymat_A_hat_times_z  DUET_NAMESPACE(polymat_A_hat_times_z)
void polymat_A_hat_times_z(polyveck *w,
                            const polyvecL A_hat[DUET_K],
                            const polyvecL *z);

#define polymat_A_hat_times_z_ntt  DUET_NAMESPACE(polymat_A_hat_times_z_ntt)
void polymat_A_hat_times_z_ntt(polyveck *w,
                                const poly  A0[DUET_K][DUET_ELL],
                                const polyveck *b,
                                const polyvecL *z);

/* ---- t = A * s1_hat  (k rows, NTT domain multiply) ---- */
#define polymat_A_times_s1_ntt  DUET_NAMESPACE(polymat_A_times_s1_ntt)
void polymat_A_times_s1_ntt(polyveck *t,
                             const poly A[DUET_K][DUET_ELL],
                             const polyvecell *s1hat);

#endif /* !DUET_POLYMAT_H */
