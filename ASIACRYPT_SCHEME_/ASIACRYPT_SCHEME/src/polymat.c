
#include "polymat.h"
#include "ntt.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "reduce.h"
#include <string.h>
#include <stdint.h>

/* Expand A0 in Z_q^{k x ell} from 32-byte seed. */
void polymat_expand_A(poly A0[DUET_K][DUET_ELL],
                      const uint8_t rho[DUET_SEEDBYTES]) {
    for (unsigned i = 0; i < DUET_K; ++i)
        for (unsigned j = 0; j < DUET_ELL; ++j)
            poly_uniform(&A0[i][j], rho, (uint16_t)((i << 8) | j));
}

/* Build A_hat in Z_{2q}^{k x L}  from A0 and b.
 *   Row i col 0      : -2*b_i + q*j  (constant term gets +q)
 *   Row i cols 1..ell: 2*A0[i][j-1]
 *   Row i col ell+1+i: 2  (identity block)
 *   All other cols   : 0
 */
void polymat_build_A_hat(polyvecL A_hat[DUET_K],
                         const poly A0[DUET_K][DUET_ELL],
                         const polyveck *b) {
    for (unsigned i = 0; i < DUET_K; ++i) {
        memset(&A_hat[i], 0, sizeof(polyvecL));

        /* Col 0: -2*b_i + q*j */
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t neg2b = (int32_t)((-2 * (int64_t)b->vec[i].coeffs[r]
                             % DUET_DQ + DUET_DQ) % DUET_DQ);
            A_hat[i].vec[0].coeffs[r] = neg2b;
        }
        /* +q on constant term for the j=1 unit vector */
        A_hat[i].vec[0].coeffs[0] =
            (int32_t)(((int64_t)A_hat[i].vec[0].coeffs[0] + DUET_Q) % DUET_DQ);

        /* Cols 1..ell: 2*A0[i][j-1] mod 2q */
        for (unsigned j = 0; j < DUET_ELL; ++j)
            for (unsigned r = 0; r < DUET_N; ++r) {
                int32_t v = 2 * A0[i][j].coeffs[r];
                if (v >= (int32_t)DUET_DQ) v -= (int32_t)DUET_DQ;
                A_hat[i].vec[1+j].coeffs[r] = v;
            }

        /* Col ell+1+i: identity block (coeff 0 = 2) */
        A_hat[i].vec[1+DUET_ELL+i].coeffs[0] = 2;
    }
}

/* w = A_hat * z  mod 2q   */
void polymat_A_hat_times_z(polyveck *w,
                            const polyvecL A_hat[DUET_K],
                            const polyvecL *z) {
    int64_t dq = (int64_t)DUET_DQ;
    for (unsigned i = 0; i < DUET_K; ++i) {
        int64_t acc[DUET_N] = {0};
        for (unsigned j = 0; j < DUET_L; ++j) {
            for (unsigned p = 0; p < DUET_N; ++p) {
                int64_t ap = A_hat[i].vec[j].coeffs[p];
                if (ap < 0) ap += dq;
                if (!ap) continue;
                for (unsigned q2 = 0; q2 < DUET_N; ++q2) {
                    int64_t bq = z->vec[j].coeffs[q2];
                    if (bq < 0) bq += dq;
                    unsigned idx = p + q2;
                    if (idx < DUET_N) acc[idx]        += ap * bq;
                    else              acc[idx-DUET_N]  -= ap * bq;
                }
            }
        }
        for (unsigned r = 0; r < DUET_N; ++r) {
            int64_t v = acc[r] % dq;
            if (v < 0) v += dq;
            w->vec[i].coeffs[r] = (int32_t)v;
        }
    }
}

void polymat_A_hat_times_z_ntt(polyveck *w,
                                const poly  A0[DUET_K][DUET_ELL],
                                const polyveck *b,
                                const polyvecL *z) {
    poly z0_ntt = z->vec[0];
    for (unsigned r = 0; r < DUET_N; ++r) {
        int32_t v = z0_ntt.coeffs[r];
        if (v < 0) v += (int32_t)DUET_Q;
        z0_ntt.coeffs[r] = v;
    }
    poly_ntt(&z0_ntt);

    /* ---- Pre-compute NTT of z[1..ell] (shared across all K rows) ---- */
    poly z_ell_ntt[DUET_ELL];
    for (unsigned j = 0; j < DUET_ELL; ++j) {
        z_ell_ntt[j] = z->vec[1 + j];
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t v = z_ell_ntt[j].coeffs[r];
            if (v < 0) v += (int32_t)DUET_Q;
            z_ell_ntt[j].coeffs[r] = v;
        }
        poly_ntt(&z_ell_ntt[j]);
    }

    /* ---- Per-row computation ---- */
    for (unsigned i = 0; i < DUET_K; ++i) {

        /* === Term 1: -2*(b_i * z[0]) mod 2q ===
         * b coeffs are in [0,q) from KeyGen.                         */
        poly b_ntt = b->vec[i];
        poly_ntt(&b_ntt);

        poly prod_b_z0;
        poly_pointwise_montgomery(&prod_b_z0, &b_ntt, &z0_ntt);
        poly_invntt_tomont(&prod_b_z0);
        poly_reduce(&prod_b_z0);   /* -> (-q, q) */

        poly neg2_b_z0;
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t v = -2 * prod_b_z0.coeffs[r];
            /* bring into [0, 2q) */
            v = (int32_t)(((int64_t)v % (int32_t)DUET_DQ
                          + (int32_t)DUET_DQ) % (int32_t)DUET_DQ);
            neg2_b_z0.coeffs[r] = v;
        }

        /* === Term 2: q * z[0]  (j*z[0] = z[0], scalar multiply by q) ===
         * z[0] coeffs are in [-B,B].  Multiply each by q, reduce mod 2q. */
        poly q_z0;
        for (unsigned r = 0; r < DUET_N; ++r) {
            int64_t v = (int64_t)DUET_Q * z->vec[0].coeffs[r];
            int32_t v2 = (int32_t)(v % (int32_t)DUET_DQ);
            if (v2 < 0) v2 += (int32_t)DUET_DQ;
            q_z0.coeffs[r] = v2;
        }

        /* === Term 3: 2 * sum_{j=0..ell-1} A0[i][j]*z[j+1] mod 2q ===
         * A0 coeffs in [0,q); z NTTs pre-computed above.             */
        poly sum_A_z;
        memset(&sum_A_z, 0, sizeof(poly));

        for (unsigned j = 0; j < DUET_ELL; ++j) {
            poly A_ntt = A0[i][j];   /* [0,q) standard domain */
            poly_ntt(&A_ntt);

            poly prod_A_z;
            poly_pointwise_montgomery(&prod_A_z, &A_ntt, &z_ell_ntt[j]);
            poly_invntt_tomont(&prod_A_z);
            poly_reduce(&prod_A_z); /* -> (-q, q) */

            poly_add(&sum_A_z, &sum_A_z, &prod_A_z);
        }

        poly two_sum_A_z;
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t v = 2 * sum_A_z.coeffs[r];
            v = (int32_t)(((int64_t)v % (int32_t)DUET_DQ
                          + (int32_t)DUET_DQ) % (int32_t)DUET_DQ);
            two_sum_A_z.coeffs[r] = v;
        }

        /* === Term 4: 2 * z[ell+1+i]  (identity-block column) ===
         * Caller zeroed these slots when z2 is absent (verify path). */
        poly two_z_id;
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t v = 2 * z->vec[1 + DUET_ELL + i].coeffs[r];
            v = (int32_t)(((int64_t)v % (int32_t)DUET_DQ
                          + (int32_t)DUET_DQ) % (int32_t)DUET_DQ);
            two_z_id.coeffs[r] = v;
        }

        /* === Accumulate all four terms mod 2q === */
        for (unsigned r = 0; r < DUET_N; ++r) {
            int64_t acc = (int64_t)neg2_b_z0.coeffs[r]
                        + (int64_t)q_z0.coeffs[r]
                        + (int64_t)two_sum_A_z.coeffs[r]
                        + (int64_t)two_z_id.coeffs[r];
            int32_t res = (int32_t)(acc % (int32_t)DUET_DQ);
            if (res < 0) res += (int32_t)DUET_DQ;
            w->vec[i].coeffs[r] = res;
        }
    }
}

/* t = A0 * s1hat  (NTT-domain pointwise multiply + iNTT)
 * s1hat must be in NTT domain (already converted by caller). */
void polymat_A_times_s1_ntt(polyveck *t,
                             const poly A0[DUET_K][DUET_ELL],
                             const polyvecell *s1hat) {
    poly A_ntt, tmp, acc;
    for (unsigned i = 0; i < DUET_K; ++i) {
        memset(&acc, 0, sizeof(poly));
        for (unsigned j = 0; j < DUET_ELL; ++j) {
            A_ntt = A0[i][j];
            poly_ntt(&A_ntt);
            poly_pointwise_montgomery(&tmp, &A_ntt, &s1hat->vec[j]);
            poly_add(&acc, &acc, &tmp);
        }
        poly_invntt_tomont(&acc);
        poly_reduce(&acc);
        t->vec[i] = acc;
    }
}
