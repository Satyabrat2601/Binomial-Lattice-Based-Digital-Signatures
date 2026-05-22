
#include "api.h"
#include "encoding.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"
#include <stdint.h>
#include <string.h>

uint32_t duet_sign_last_attempts = 0;

/* ================================================================
   Helpers
   ================================================================ */

/* Build S = (j | s1 | s2)^T as a polyvecL.
 * S[0] = j = constant polynomial 1.
 * S[1..ell] = s1[0..ell-1].
 * S[ell+1..L-1] = s2[0..k-1]. */
static void _build_S(polyvecL *S, const polyvecell *s1, const polyveck *s2) {
    memset(&S->vec[0], 0, sizeof(poly));
    S->vec[0].coeffs[0] = 1;
    for (unsigned i = 0; i < DUET_ELL; ++i) S->vec[1+i]          = s1->vec[i];
    for (unsigned i = 0; i < DUET_K;   ++i) S->vec[1+DUET_ELL+i] = s2->vec[i];
}

/* Lift signed z-vector into Z_{2q}: add 2q to negative coefficients. */
static void _lift_to_2q(polyvecL *out, const polyvecL *in) {
    for (unsigned i = 0; i < DUET_L; ++i)
        for (unsigned r = 0; r < DUET_N; ++r) {
            int32_t c = in->vec[i].coeffs[r];
            out->vec[i].coeffs[r] = c < 0 ? c + (int32_t)DUET_DQ : c;
        }
}

static void _expand_challenge_from_seed(poly *c, const uint8_t seed[DUET_SEEDBYTES]) {
    keccak_state cst;
    uint8_t pos[DUET_N];    /* position bag: pos[i] = actual position i represents */
    uint8_t outbuf[DUET_N]; /* 1 draw byte per position */

    for (unsigned i = 0; i < DUET_N; ++i) pos[i] = (uint8_t)i;

    xof256_absorb_once(&cst, seed, DUET_SEEDBYTES);
    xof256_squeeze(outbuf, DUET_N, &cst);

    memset(c->coeffs, 0, sizeof(c->coeffs));
    for (unsigned k = 0; k < DUET_TAU; ++k) {
        unsigned range = DUET_N - k;
        unsigned j = (unsigned)outbuf[k] % range;  /* uniform mod N (power of 2, unbiased) */
        c->coeffs[pos[j]] = 1;
        pos[j] = pos[range - 1]; /* remove from bag */
    }
}

static void _pack_w1_for_hash(uint8_t *buf, size_t *len_out, const polyveck *w) {
    polyveck_pack_hi(buf, w);
    *len_out = DUET_K * DUET_N;
}

/* ================================================================
   KeyGen  (Fig.1)
   ================================================================ */
int crypto_sign_keypair_internal(uint8_t *pk, uint8_t *sk,
                                  const uint8_t seed[DUET_SEEDBYTES]) {
    uint8_t seedbuf[2*DUET_SEEDBYTES + DUET_CRHBYTES];
    uint8_t rho[DUET_SEEDBYTES], sigma[DUET_CRHBYTES], key[DUET_SEEDBYTES];
    uint8_t tr[DUET_SEEDBYTES];
    keccak_state st;

    xof256_absorb_once(&st, seed, DUET_SEEDBYTES);
    xof256_squeeze(seedbuf, sizeof(seedbuf), &st);
    memcpy(rho,   seedbuf,                            DUET_SEEDBYTES);
    memcpy(sigma, seedbuf + DUET_SEEDBYTES,           DUET_CRHBYTES);
    memcpy(key,   seedbuf + DUET_SEEDBYTES + DUET_CRHBYTES, DUET_SEEDBYTES);

    /* Step 1: A0 <- R_q^{k x ell} */
    poly A0[DUET_K][DUET_ELL];
    polymat_expand_A(A0, rho);

    /* Step 2: (s1, s2) <- CBD_eta */
    polyvecell s1; polyveck s2;
    polyvecell_expand_eta(&s1, sigma, 0);
    polyveck_expand_eta(&s2, sigma, (uint16_t)DUET_ELL);

    /* Step 3: b = A0*s1 + s2 mod q (NTT multiplication) */
    polyvecell s1h = s1; polyvecell_ntt(&s1h);
    polyveck b;
    polymat_A_times_s1_ntt(&b, A0, &s1h);
    polyveck_add(&b, &b, &s2); polyveck_reduce2q(&b); polyveck_freeze(&b);

    pack_pk(pk, rho, &b);
    xof256_absorb_once(&st, pk, DUET_CRYPTO_PUBLICKEYBYTES);
    xof256_squeeze(tr, DUET_SEEDBYTES, &st);
    pack_sk(sk, rho, key, tr, &s1, &s2);
    return 0;
}

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
    uint8_t seed[DUET_SEEDBYTES];
    if (randombytes(seed, DUET_SEEDBYTES) != 0) return -1;
    return crypto_sign_keypair_internal(pk, sk, seed);
}

/* ================================================================
   Sign core  (Fig.2 CSign rejection loop)
   Returns: z1_out, h_out, c_seed_out on success.
   ================================================================ */
static int _sign_core(
    polyvecL       *z1_out,
    polyveck       *h_out,
    uint8_t         c_seed_out[DUET_SEEDBYTES],
    uint32_t       *attempts_out,
    const polyvecL  A_hat[DUET_K],
    const polyvecL *S,
    const uint8_t  *mu,
    const uint8_t  *mask_seed)
{
    polyvecL y, y_2q, z, z1;
    polyveck z2, w, wp;
    poly c;
    uint8_t w1buf[DUET_K * DUET_N];
    uint8_t hash_in[DUET_K * DUET_N + DUET_SEEDBYTES];
    size_t  w1len;
    uint8_t y_seed[DUET_SEEDBYTES];
    uint8_t ab[4];
    keccak_state st;
    uint32_t attempt = 0;

    for (;; ++attempt) {
        /* Per-attempt y seed: derived deterministically from mask_seed + attempt counter */
        ab[0]=(uint8_t)attempt; ab[1]=(uint8_t)(attempt>>8);
        ab[2]=(uint8_t)(attempt>>16); ab[3]=(uint8_t)(attempt>>24);
        xof256_absorb_twice(&st, mask_seed, DUET_SEEDBYTES, ab, 4);
        xof256_squeeze(y_seed, DUET_SEEDBYTES, &st);

        /* Fig.1 step 1: y <- CBD_{gamma1}^L */
        polyvecL_expand_cbd_gamma1(&y, y_seed);

        /* Fig.1 step 2: w = A_hat * y mod 2q */
        _lift_to_2q(&y_2q, &y);
        polymat_A_hat_times_z(&w, A_hat, &y_2q);

        _pack_w1_for_hash(w1buf, &w1len, &w);
        memcpy(hash_in, w1buf, w1len);
        memcpy(hash_in + w1len, mu, DUET_SEEDBYTES);
        xof256_absorb_once(&st, hash_in, w1len + DUET_SEEDBYTES);
        xof256_squeeze(c_seed_out, DUET_SEEDBYTES, &st);

        _expand_challenge_from_seed(&c, c_seed_out);

        polyvecL_mul_xiS(&z, &y, &c, S);

        /* Rejection: ||z||_inf > B */
        if (polyvecL_chknorm(&z, (int32_t)DUET_B)) continue;

        /* Rejection: z not in CBD_B (statistical validity) */
        if (!polyvecL_is_in_cbd_bound(&z, (int32_t)DUET_CBD_BOUND(DUET_B))) continue;

        /* Fig.2 step 9: Split(z) -> (z1, z2) */
        polyvecL_split(&z1, &z2, &z);

        /* Fig.2 step 10: wp = w - 2*z2 mod 2q */
        polyveck_sub_2z2(&wp, &w, &z2);

        /* Fig.2 step 10: h = MakeHint(w, wp) = (HighBits(w) - HighBits(wp)) mod H_RANGE */
        polyveck_make_hint(h_out, &w, &wp);

        *z1_out = z1;
        if (attempts_out) *attempts_out = attempt + 1;
        duet_sign_last_attempts = attempt + 1;
        return 0;
    }
}

/* ================================================================
   Sign  (Fig.1 / Fig.2)
   ================================================================ */
int crypto_sign_signature_internal(
    uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *pre, size_t prelen,
    const uint8_t rnd[DUET_SEEDBYTES],
    const uint8_t *sk)
{
    uint8_t rho[DUET_SEEDBYTES], key[DUET_SEEDBYTES], tr[DUET_SEEDBYTES];
    uint8_t mu[DUET_SEEDBYTES], mask_seed[DUET_SEEDBYTES];
    uint8_t c_seed[DUET_SEEDBYTES];
    poly A0[DUET_K][DUET_ELL];
    polyvecell s1; polyveck s2, b;
    polyvecL A_hat[DUET_K], S, z1;
    polyveck h;
    keccak_state st;

    unpack_sk(rho, key, tr, &s1, &s2, sk);

    polymat_expand_A(A0, rho);
    { polyvecell s1h = s1; polyvecell_ntt(&s1h);
      polymat_A_times_s1_ntt(&b, A0, &s1h); }
    polyveck_add(&b, &b, &s2); polyveck_reduce2q(&b); polyveck_freeze(&b);
    polymat_build_A_hat(A_hat, A0, &b);
    _build_S(&S, &s1, &s2);

    /* mu = H(tr || pre || m) */
    xof256_absorb_thrice(&st, tr, DUET_SEEDBYTES, pre, prelen, m, mlen);
    xof256_squeeze(mu, DUET_SEEDBYTES, &st);

    /* mask_seed = H(key || rnd || mu) */
    xof256_absorb_thrice(&st, key, DUET_SEEDBYTES, rnd, DUET_SEEDBYTES, mu, DUET_SEEDBYTES);
    xof256_squeeze(mask_seed, DUET_SEEDBYTES, &st);

    uint32_t attempts;
    _sign_core(&z1, &h, c_seed, &attempts, A_hat, &S, mu, mask_seed);

    pack_sig(sig, &z1, &h, c_seed);
    *siglen = DUET_CRYPTO_BYTES;
    return 0;
}

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen,
    const uint8_t *sk)
{
    if (ctxlen > 255) return -1;
    uint8_t pre[256]; pre[0] = (uint8_t)ctxlen;
    if (ctxlen && ctx) memcpy(pre+1, ctx, ctxlen);
    uint8_t rnd[DUET_SEEDBYTES];
    if (randombytes(rnd, DUET_SEEDBYTES) != 0) return -1;
    return crypto_sign_signature_internal(sig, siglen, m, mlen,
                                          pre, 1+ctxlen, rnd, sk);
}

int crypto_sign(uint8_t *sm, size_t *smlen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk)
{
    uint8_t sig[DUET_CRYPTO_BYTES];   
    size_t  siglen;
    int r = crypto_sign_signature(sig, &siglen, m, mlen, ctx, ctxlen, sk);
    if (r) return r;
    /* Now m is untouched; move it to make room for the signature */
    memmove(sm + siglen, m, mlen);         
    memcpy(sm, sig, siglen);                
    *smlen = siglen + mlen;
    return 0;
}

int crypto_sign_verify_internal(
    const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *pre, size_t prelen,
    const uint8_t *pk)
{
    if (siglen != DUET_CRYPTO_BYTES) return -1;

    uint8_t rho[DUET_SEEDBYTES], tr[DUET_SEEDBYTES], mu[DUET_SEEDBYTES];
    uint8_t c_seed[DUET_SEEDBYTES], c_seed_prime[DUET_SEEDBYTES];
    uint8_t w1buf[DUET_K * DUET_N];
    uint8_t hash_in[DUET_K * DUET_N + DUET_SEEDBYTES];
    size_t  w1len;
    poly A0[DUET_K][DUET_ELL];
    polyveck b, w_prime, h, hi_rec;
    polyvecL z1;
    keccak_state st;

    unpack_pk(rho, &b, pk);
    polymat_expand_A(A0, rho);

    /* Zero-initialise z1 fully so the z2 block (vec[ell+1..L-1]),
     * which is NOT transmitted in the signature, stays zero.
     * polymat_A_hat_times_z_ntt reads those slots for the identity column. */
    memset(&z1, 0, sizeof(z1));

    /* Unpack sig -> (z1[0..ell], h, c_seed) */
    if (unpack_sig(&z1, &h, c_seed, sig) != 0) return -1;

    /* Norm check on z1 (ell+1 polys, signed coefficients in [-B,B]) */
    for (unsigned i = 0; i <= DUET_ELL; ++i)
        if (poly_chknorm(&z1.vec[i], (int32_t)DUET_B)) return -1;

    /* Re-derive tr = H(pk) and mu = H(tr || pre || m) */
    xof256_absorb_once(&st, pk, DUET_CRYPTO_PUBLICKEYBYTES);
    xof256_squeeze(tr, DUET_SEEDBYTES, &st);
    xof256_absorb_thrice(&st, tr, DUET_SEEDBYTES, pre, prelen, m, mlen);
    xof256_squeeze(mu, DUET_SEEDBYTES, &st);

    poly c;
    _expand_challenge_from_seed(&c, c_seed);

    polymat_A_hat_times_z_ntt(&w_prime, A0, &b, &z1);
    for (unsigned r = 0; r < DUET_N; ++r) {
        if (c.coeffs[r] == 0) continue;
        for (unsigned i = 0; i < DUET_K; ++i) {
            int32_t v = w_prime.vec[i].coeffs[r] - (int32_t)DUET_Q;
            if (v < 0) v += (int32_t)DUET_DQ;
            w_prime.vec[i].coeffs[r] = v;
        }
    }

    for (unsigned i = 0; i < DUET_K; ++i)
        poly_usehint(&hi_rec.vec[i], &w_prime.vec[i], &h.vec[i]);

    for (unsigned i = 0; i < DUET_K; ++i)
        for (unsigned r = 0; r < DUET_N; ++r)
            w1buf[i * DUET_N + r] = (uint8_t)(hi_rec.vec[i].coeffs[r] & 0xFF);
    w1len = DUET_K * DUET_N;

    memcpy(hash_in, w1buf, w1len);
    memcpy(hash_in + w1len, mu, DUET_SEEDBYTES);
    xof256_absorb_once(&st, hash_in, w1len + DUET_SEEDBYTES);
    xof256_squeeze(c_seed_prime, DUET_SEEDBYTES, &st);

    /* Constant-time comparison */
    uint8_t diff = 0;
    for (unsigned i = 0; i < DUET_SEEDBYTES; ++i)
        diff |= c_seed[i] ^ c_seed_prime[i];
    return diff ? -1 : 0;
}

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk)
{
    if (ctxlen > 255) return -1;
    uint8_t pre[256]; pre[0] = (uint8_t)ctxlen;
    if (ctxlen && ctx) memcpy(pre+1, ctx, ctxlen);
    return crypto_sign_verify_internal(sig, siglen, m, mlen, pre, 1+ctxlen, pk);
}

int crypto_sign_open(uint8_t *m, size_t *mlen,
    const uint8_t *sm, size_t smlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk)
{
    if (smlen < DUET_CRYPTO_BYTES) { *mlen = 0; return -1; }
    *mlen = smlen - DUET_CRYPTO_BYTES;
    if (crypto_sign_verify(sm, DUET_CRYPTO_BYTES,
                           sm + DUET_CRYPTO_BYTES, *mlen,
                           ctx, ctxlen, pk) != 0) { *mlen = 0; return -1; }
    memmove(m, sm + DUET_CRYPTO_BYTES, *mlen);
    return 0;
}

/* ================================================================
   CSign / CVerify  (Fig.2 -- same core)
   ================================================================ */
int crypto_sign_signature_compact(uint8_t *sig, size_t *siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk)
{ return crypto_sign_signature(sig, siglen, m, mlen, ctx, ctxlen, sk); }

int crypto_sign_verify_compact(const uint8_t *sig, size_t siglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk)
{ return crypto_sign_verify(sig, siglen, m, mlen, ctx, ctxlen, pk); }

int crypto_sign_signature_rans(uint8_t *rsig, size_t *rsiglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk)
{
    uint8_t sig[DUET_CRYPTO_BYTES]; size_t siglen;
    if (crypto_sign_signature(sig, &siglen, m, mlen, ctx, ctxlen, sk)) return -1;
    encode_signature_rans(rsig, rsiglen, sig, siglen);
    return 0;
}

int crypto_sign_verify_rans(const uint8_t *rsig, size_t rsiglen,
    const uint8_t *m, size_t mlen,
    const uint8_t *ctx, size_t ctxlen, const uint8_t *pk)
{
    uint8_t sig[DUET_CRYPTO_BYTES]; size_t siglen;
    if (decode_signature_rans(sig, &siglen, rsig, rsiglen)) return -1;
    return crypto_sign_verify(sig, siglen, m, mlen, ctx, ctxlen, pk);
}
