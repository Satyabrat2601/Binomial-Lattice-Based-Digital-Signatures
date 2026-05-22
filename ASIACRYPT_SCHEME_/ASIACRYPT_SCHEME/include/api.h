
#ifndef DUET_API_H
#define DUET_API_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

#define CRYPTO_ALGNAME              DUET_CRYPTO_ALGNAME
#define CRYPTO_PUBLICKEYBYTES       DUET_CRYPTO_PUBLICKEYBYTES
#define CRYPTO_SECRETKEYBYTES       DUET_CRYPTO_SECRETKEYBYTES
#define CRYPTO_BYTES                DUET_CRYPTO_BYTES
#define CRYPTO_CSIG_BYTES           DUET_CRYPTO_BYTES  /* compact sig = same format */
/* rANS sig: z1_rans + h_verbatim + c_seed. z1 compresses from 10 bits to ~5.75 bits/coeff. */
#define CRYPTO_RANS_SIG_MAXBYTES    DUET_CRYPTO_BYTES  /* safe upper bound */

/* NIST PQC API */
#define crypto_sign_keypair          DUET_NAMESPACE(keypair)
#define crypto_sign_signature        DUET_NAMESPACE(signature)
#define crypto_sign                  DUET_NAMESPACE(sign)
#define crypto_sign_verify           DUET_NAMESPACE(verify)
#define crypto_sign_open             DUET_NAMESPACE(open)

/* Internal (for testing / KATs) */
#define crypto_sign_keypair_internal       DUET_NAMESPACE(keypair_internal)
#define crypto_sign_signature_internal     DUET_NAMESPACE(signature_internal)
#define crypto_sign_verify_internal        DUET_NAMESPACE(verify_internal)
#define crypto_sign_signature_compact      DUET_NAMESPACE(signature_compact)
#define crypto_sign_verify_compact         DUET_NAMESPACE(verify_compact)
#define crypto_sign_signature_rans         DUET_NAMESPACE(signature_rans)
#define crypto_sign_verify_rans            DUET_NAMESPACE(verify_rans)

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);
int crypto_sign_keypair_internal(uint8_t *pk, uint8_t *sk,
                                  const uint8_t seed[DUET_SEEDBYTES]);

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen,
    const uint8_t *sk);
int crypto_sign_signature_internal(uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *pre, size_t prelen,
    const uint8_t rnd[DUET_SEEDBYTES],
    const uint8_t *sk);

int crypto_sign(uint8_t *sm, size_t *smlen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen,
    const uint8_t *sk);

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen,
    const uint8_t *pk);
int crypto_sign_verify_internal(const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *pre, size_t prelen,
    const uint8_t *pk);

int crypto_sign_open(uint8_t *m, size_t *mlen,
    const uint8_t *sm, size_t smlen,
    const uint8_t *ctx, size_t ctxlen,
    const uint8_t *pk);

int crypto_sign_signature_compact(uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk);
int crypto_sign_verify_compact(const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk);
int crypto_sign_signature_rans(uint8_t *rsig, size_t *rsiglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk);
int crypto_sign_verify_rans(const uint8_t *rsig, size_t rsiglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk);

/* Attempt counter -- set by each sign call */
extern uint32_t duet_sign_last_attempts;

#endif /* DUET_API_H */
