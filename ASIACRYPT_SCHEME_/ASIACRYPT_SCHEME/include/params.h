
#ifndef DUET_PARAMS_H
#define DUET_PARAMS_H

#include "config.h"
#include <stdint.h>

/* Namespace prefix to allow multiple modes in one binary */
#define DUET_NAMESPACE(s) cryptolab_duetsign_##s

/* ============================================================ Shared */
#define DUET_N            256
#define DUET_Q          64513
#define DUET_DQ        129026   /* 2*q */
#define DUET_SEEDBYTES     32
#define DUET_CRHBYTES      64
#define DUET_ETA            1
#define DUET_B_BITS        10   /* bits for z packing: 2^10=1024 >= 2*450+1=901 */

/* ============================================================
 * DUET-I  (mode1): ell=3 k=2 eta=1 B=300 tau=39 omega=5
 * Security ~120-bit
 * Sig: 4*320 + 2*256 + 32 = 1280+512+32 = 1824 B
 * ============================================================ */
#if DUET_CONFIG_MODE == DUET_MODE1

#define DUET_K              2
#define DUET_ELL            3
#define DUET_TAU           39
#define DUET_BETA          39   /* tau*eta */
#define DUET_B            300
#define DUET_GAMMA1       339   /* B+beta */
#define DUET_ALPHA        678   /* 2*gamma1 */
#define DUET_OMEGA          5
#define DUET_L              (DUET_K + DUET_ELL + 1)  /* 6 */
/* HighBits index range = ceil(2q / alpha) = ceil(190.3) = 191 */
#define DUET_H_RANGE      191
#define DUET_CRYPTO_ALGNAME "duet-sign-mode1"

/* ============================================================
 * DUET-II (mode3): ell=4 k=3 eta=1 B=420 tau=49 omega=5
 * Security ~180-bit 
 * Sig: 5*320 + 3*256 + 32 = 1600+768+32 = 2400 B
 * ============================================================ */
#elif DUET_CONFIG_MODE == DUET_MODE3

#define DUET_K              3
#define DUET_ELL            4
#define DUET_TAU           49
#define DUET_BETA          49
#define DUET_B            420
#define DUET_GAMMA1       469
#define DUET_ALPHA        938
#define DUET_OMEGA          5
#define DUET_L              (DUET_K + DUET_ELL + 1)  /* 8 */
/* ceil(129026/938) = ceil(137.6) = 138 */
#define DUET_H_RANGE      138
#define DUET_CRYPTO_ALGNAME "duet-sign-mode3"

/* ============================================================
 * DUET-III(mode5): ell=6 k=4 eta=1 B=450 tau=60 omega=5
 * Security ~260-bit 
 * Sig: 7*320 + 4*256 + 32 = 2240+1024+32 = 3296 B
 * ============================================================ */
#elif DUET_CONFIG_MODE == DUET_MODE5

#define DUET_K              4
#define DUET_ELL            6
#define DUET_TAU           60
#define DUET_BETA          60
#define DUET_B            450
#define DUET_GAMMA1       510
#define DUET_ALPHA       1020
#define DUET_OMEGA          5
#define DUET_L              (DUET_K + DUET_ELL + 1)  /* 11 */
/* ceil(129026/1020) = ceil(126.5) = 127 */
#define DUET_H_RANGE      127
#define DUET_CRYPTO_ALGNAME "duet-sign-mode5"

#else
#error "Unknown DUET_CONFIG_MODE"
#endif

/* ============================================================ Derived sizes */
/* CBD bound: 85% of B (used in rejection sampling) */
#define DUET_CBD_BOUND(B)    (((B) * 85) / 100)

/* Polynomial packing sizes */
#define DUET_POLYQ_PACKEDBYTES    512   /* 16 bits/coeff * 256 = 512 B */
#define DUET_POLYETA_PACKEDBYTES   64   /* 2 bits/coeff * 256 / 8 = 64 B */
#define DUET_POLYZ_PACKEDBYTES    320   /* B_BITS=10 bits/coeff * 256 / 8 = 320 B */
#define DUET_POLYH_PACKEDBYTES    256   /* 8 bits/coeff * 256 = 256 B (h in [0,H_RANGE)) */
#define DUET_POLYC_PACKEDBYTES     32   /* c_seed = 32-byte SHAKE hash */

/* Key sizes (match estimator exactly) */
#define DUET_CRYPTO_PUBLICKEYBYTES \
    (DUET_SEEDBYTES + DUET_K * DUET_POLYQ_PACKEDBYTES)

#define DUET_CRYPTO_SECRETKEYBYTES \
    (3*DUET_SEEDBYTES + (DUET_ELL + DUET_K) * DUET_POLYETA_PACKEDBYTES)

/* Signature size: z1 || h || c_seed */
#define DUET_CRYPTO_BYTES \
    ((DUET_ELL + 1) * DUET_POLYZ_PACKEDBYTES \
     + DUET_K * DUET_POLYH_PACKEDBYTES \
     + DUET_POLYC_PACKEDBYTES)

/* Sanity checks */
#if (DUET_GAMMA1 != (DUET_B + DUET_BETA))
#error "gamma1 must equal B + beta"
#endif
#if (DUET_ALPHA != (2 * DUET_GAMMA1))
#error "alpha must equal 2*gamma1"
#endif

#endif /* DUET_PARAMS_H */
