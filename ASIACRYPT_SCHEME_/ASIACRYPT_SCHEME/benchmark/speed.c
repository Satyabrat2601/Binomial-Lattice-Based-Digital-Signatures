
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "api.h"
#include "cpucycles.h"
#include "params.h"
#include "randombytes.h"
#include "speed_print.h"

#define NTESTS  1000
#define MLEN    DUET_SEEDBYTES

static uint64_t t[NTESTS + 1];

int main(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t csig[CRYPTO_CSIG_BYTES];
    uint8_t rsig[CRYPTO_RANS_SIG_MAXBYTES];
    uint8_t msg[MLEN];
    size_t  siglen, csiglen, rsiglen;
    clock_t cts0, cts1;
    int i;

    printf("\n");
    printf("========================================================\n");
    printf("  DUET-Sign Benchmark  --  %s\n", CRYPTO_ALGNAME);
    printf("========================================================\n");
    printf("  Parameters:\n");
    printf("    n=%-3d  k=%-d  ell=%-d  q=%-6d  tau=%-3d  B=%-6d\n",
           DUET_N, DUET_K, DUET_ELL, DUET_Q, DUET_TAU, DUET_B);
    printf("    gamma1=%-6d  eta=%d  H_RANGE=%d  H_BITS=%d\n",
           DUET_GAMMA1, DUET_ETA, DUET_H_RANGE, DUET_H_BITS);
    printf("  Key sizes:\n");
    printf("    pk=%d B   sk=%d B\n",
           CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES);
    printf("  Signature sizes:\n");
    printf("    Plain:      %d B\n", CRYPTO_BYTES);
    printf("    CSign:      %d B  (%.1f%% of plain)\n",
           CRYPTO_CSIG_BYTES, 100.0*CRYPTO_CSIG_BYTES/CRYPTO_BYTES);
    printf("    rANS max:   %d B  (variable-length)\n", CRYPTO_RANS_SIG_MAXBYTES);
    printf("  NTESTS=%d  MLEN=%d\n", NTESTS, MLEN);
    printf("========================================================\n\n");
    printf("  %-44s  %8s  %8s  %s\n", "Operation", "median", "avg", "cycles");
    printf("  %s\n", "------------------------------------------------------------");

    randombytes(msg, MLEN);

    /* ---- KeyGen ---- */
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_keypair(pk, sk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_keypair", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS);
    printf("\n");

    /* ---- Plain Sign ---- */
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_signature(sig, &siglen, msg, MLEN, NULL, 0, sk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_signature (plain)", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op  sig=%zu B\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS, siglen);
    printf("\n");

    /* ---- Plain Verify ---- */
    crypto_sign_signature(sig, &siglen, msg, MLEN, NULL, 0, sk);
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_verify(sig, siglen, msg, MLEN, NULL, 0, pk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_verify   (plain)", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS);
    printf("\n");

    /* ---- CSign Sign ---- */
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_signature_compact(csig, &csiglen, msg, MLEN, NULL, 0, sk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_signature_compact", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op  csig=%zu B\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS, csiglen);
    printf("\n");

    /* ---- CSign Verify ---- */
    crypto_sign_signature_compact(csig, &csiglen, msg, MLEN, NULL, 0, sk);
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_verify_compact(csig, csiglen, msg, MLEN, NULL, 0, pk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_verify_compact   ", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS);
    printf("\n");

    /* ---- rANS Sign ---- */
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_signature_rans(rsig, &rsiglen, msg, MLEN, NULL, 0, sk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_signature_rans   ", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op  rsig=%zu B\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS, rsiglen);
    printf("\n");

    /* ---- rANS Verify ---- */
    crypto_sign_signature_rans(rsig, &rsiglen, msg, MLEN, NULL, 0, sk);
    cts0 = clock();
    for (i = 0; i <= NTESTS; ++i) {
        t[i] = cpucycles();
        crypto_sign_verify_rans(rsig, rsiglen, msg, MLEN, NULL, 0, pk);
    }
    cts1 = clock();
    printf("  "); print_results("crypto_sign_verify_rans      ", t, NTESTS+1);
    printf("  %44s  wall: %.3f ms/op\n", "",
           (double)(cts1-cts0)*1000.0/CLOCKS_PER_SEC/NTESTS);
    printf("\n");

    /* ---- Summary ---- */
    printf("========================================================\n");
    printf("  Signature size comparison:\n");
    printf("    Plain  : %5d B   (baseline)\n", CRYPTO_BYTES);
    printf("    CSign  : %5d B   (%.1f%% smaller)\n",
           CRYPTO_CSIG_BYTES,
           100.0*(CRYPTO_BYTES - CRYPTO_CSIG_BYTES)/CRYPTO_BYTES);
    {
        crypto_sign_signature_rans(rsig, &rsiglen, msg, MLEN, NULL, 0, sk);
        /* Average over multiple messages for rANS size */
        size_t total = rsiglen;
        for (int j = 1; j < 10; ++j) {
            randombytes(msg, MLEN);
            crypto_sign_keypair(pk, sk);
            crypto_sign_signature_rans(rsig, &rsiglen, msg, MLEN, NULL, 0, sk);
            total += rsiglen;
        }
        size_t avg_rans = total / 10;
        printf("    rANS   : ~%4zu B   (%.1f%% smaller, variable)\n",
               avg_rans,
               100.0*(CRYPTO_BYTES - (int)avg_rans)/CRYPTO_BYTES);
    }
    printf("========================================================\n\n");

    return 0;
}
