

#include "api.h"
#include "params.h"
#include "encoding.h"
#include "randombytes.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#define ITERATIONS  100
#define MSG_LEN     DUET_SEEDBYTES
#define BENCH_N     10000
#define REJ_N       1000

/* ---- cycle counter ---- */
#if defined(__x86_64__) || defined(__i386__)
static inline uint64_t cpucycles(void) {
    uint32_t lo, hi;
    __asm__ volatile ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}
#else
static inline uint64_t cpucycles(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
}
#endif

static int cmp_u64(const void *a, const void *b) {
    uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b;
    return (x>y)-(x<y);
}
static void setup(uint8_t *pk, uint8_t *sk, uint8_t *msg) {
    crypto_sign_keypair(pk, sk);
    randombytes(msg, MSG_LEN);
}

/* ---- Tests ---- */
static int test_plain(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES], msg[MSG_LEN]; size_t sl;
    for (int i=0; i<ITERATIONS; ++i) {
        setup(pk, sk, msg);
        if (crypto_sign_signature(sig, &sl, msg, MSG_LEN, NULL, 0, sk)) { printf("  plain: sign fail\n"); return 1; }
        if (crypto_sign_verify(sig, sl, msg, MSG_LEN, NULL, 0, pk))      { printf("  plain: verify fail\n"); return 1; }
    }
    return 0;
}
static int test_compact(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES], msg[MSG_LEN]; size_t sl;
    for (int i=0; i<ITERATIONS; ++i) {
        setup(pk, sk, msg);
        if (crypto_sign_signature_compact(sig, &sl, msg, MSG_LEN, NULL, 0, sk)) { printf("  compact: sign fail\n"); return 1; }
        if (crypto_sign_verify_compact(sig, sl, msg, MSG_LEN, NULL, 0, pk))      { printf("  compact: verify fail\n"); return 1; }
    }
    return 0;
}
static int test_rans(void) {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t rsig[CRYPTO_RANS_SIG_MAXBYTES], msg[MSG_LEN]; size_t rsl;
    for (int i=0; i<ITERATIONS; ++i) {
        setup(pk, sk, msg);
        if (crypto_sign_signature_rans(rsig, &rsl, msg, MSG_LEN, NULL, 0, sk)) { printf("  rans: sign fail\n"); return 1; }
        if (rsl > (size_t)CRYPTO_RANS_SIG_MAXBYTES) { printf("  rans: size %zu > %d\n", rsl, CRYPTO_RANS_SIG_MAXBYTES); return 1; }
        if (crypto_sign_verify_rans(rsig, rsl, msg, MSG_LEN, NULL, 0, pk))      { printf("  rans: verify fail\n"); return 1; }
    }
    return 0;
}
static int test_cross(void) {
    uint8_t pk1[CRYPTO_PUBLICKEYBYTES], sk1[CRYPTO_SECRETKEYBYTES];
    uint8_t pk2[CRYPTO_PUBLICKEYBYTES], sk2[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES], m1[MSG_LEN], m2[MSG_LEN]; size_t sl;
    setup(pk1, sk1, m1); setup(pk2, sk2, m2);
    crypto_sign_signature(sig, &sl, m1, MSG_LEN, NULL, 0, sk1);
    if (!crypto_sign_verify(sig, sl, m1, MSG_LEN, NULL, 0, pk2)) { printf("  wrong-key not rejected\n"); return 1; }
    if (!crypto_sign_verify(sig, sl, m2, MSG_LEN, NULL, 0, pk1)) { printf("  wrong-msg not rejected\n"); return 1; }
    sig[0] ^= 1;
    if (!crypto_sign_verify(sig, sl, m1, MSG_LEN, NULL, 0, pk1)) { printf("  tampered not rejected\n"); return 1; }
    return 0;
}



/* ---- Cycle benchmark ---- */
static void benchmark_cycles(void) {
    printf("=============================================================\n");
    printf("  Cycle-Count Benchmark  (%d iterations)\n", BENCH_N);
   
    printf("=============================================================\n");

    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES], msg[MSG_LEN]; size_t sl;
    uint64_t *t=(uint64_t*)malloc(BENCH_N*sizeof(uint64_t));
    if(!t){printf("  [malloc fail]\n");return;}

    /* KeyGen */
    for(int i=0;i<BENCH_N;++i){
        uint8_t s[DUET_SEEDBYTES]; s[0]=(uint8_t)i; s[1]=(uint8_t)(i>>8);
        uint64_t c0=cpucycles(); crypto_sign_keypair_internal(pk,sk,s); t[i]=cpucycles()-c0;
    }
    qsort(t,BENCH_N,sizeof(uint64_t),cmp_u64);
    uint64_t kg_med=t[BENCH_N/2], kg_avg=0;
    for(int i=0;i<BENCH_N;i++) kg_avg+=t[i]; kg_avg/=BENCH_N;

    /* Sign */
    crypto_sign_keypair(pk,sk);
    for(int i=0;i<BENCH_N;++i){
        msg[0]=(uint8_t)i; msg[1]=(uint8_t)(i>>8);
        uint64_t c0=cpucycles(); crypto_sign_signature(sig,&sl,msg,MSG_LEN,NULL,0,sk); t[i]=cpucycles()-c0;
    }
    qsort(t,BENCH_N,sizeof(uint64_t),cmp_u64);
    uint64_t sg_med=t[BENCH_N/2], sg_avg=0;
    for(int i=0;i<BENCH_N;i++) sg_avg+=t[i]; sg_avg/=BENCH_N;

    /* Verify */
    msg[0]=(uint8_t)(BENCH_N-1); msg[1]=(uint8_t)((BENCH_N-1)>>8);
    crypto_sign_signature(sig,&sl,msg,MSG_LEN,NULL,0,sk);
    for(int i=0;i<BENCH_N;++i){
        uint64_t c0=cpucycles(); crypto_sign_verify(sig,sl,msg,MSG_LEN,NULL,0,pk); t[i]=cpucycles()-c0;
    }
    qsort(t,BENCH_N,sizeof(uint64_t),cmp_u64);
    uint64_t vf_med=t[BENCH_N/2], vf_avg=0;
    for(int i=0;i<BENCH_N;i++) vf_avg+=t[i]; vf_avg/=BENCH_N;
    free(t);

    /* Table output  */
    printf("  %-14s %12s %14s %12s\n", "Scheme", "KeyGen", "Sign", "Verify");
    printf("  %-14s %12s %14s %12s\n", "------", "------", "----", "------");
    printf("  %-10s  med %12llu %14llu %12llu\n",
           CRYPTO_ALGNAME,
           (unsigned long long)kg_med,
           (unsigned long long)sg_med,
           (unsigned long long)vf_med);
    printf("  %-10s  ave %12llu %14llu %12llu\n",
           "",
           (unsigned long long)kg_avg,
           (unsigned long long)sg_avg,
           (unsigned long long)vf_avg);
    printf("\n");
    printf("  KeyGen  : %llu cycles  (%.4f ms)\n", (unsigned long long)kg_med, kg_med/3e9*1000);
    printf("  Sign    : %llu cycles  (%.4f ms)  [incl. rejection]\n", (unsigned long long)sg_med, sg_med/3e9*1000);
    printf("  Verify  : %llu cycles  (%.4f ms)\n", (unsigned long long)vf_med, vf_med/3e9*1000);
    printf("\n");
    /* rANS size */
    {
        uint8_t rsig[CRYPTO_RANS_SIG_MAXBYTES]; size_t rsl=0;
        crypto_sign_signature_rans(rsig,&rsl,msg,MSG_LEN,NULL,0,sk);
        printf("  Sig_packed : %d B   (reference implementation)\n", CRYPTO_BYTES);
        printf("  Sig_rANS   : %zu B  (rANS entropy-coded, this implementation)\n", rsl);
        
               
               
        
        printf("  Security   : %s\n",
               (DUET_K==2)?"~120-bit ":(DUET_K==3)?"~180-bit ":"~260-bit ");
    }
    printf("=============================================================\n\n");
}

/* ---- main ---- */
int main(void) {
    printf("\n=============================================================\n");
    printf("  DUET-Sign Reference Implementation  -- ");
    printf("=============================================================\n");
    printf("  Scheme  : %s\n", CRYPTO_ALGNAME);
    printf("  n=%d  q=%d  ell=%d  k=%d  eta=%d\n", DUET_N, DUET_Q, DUET_ELL, DUET_K, DUET_ETA);
    printf("  tau=%d  B=%d  gamma1=%d  beta=%d  omega=%d\n",
           DUET_TAU, DUET_B, DUET_GAMMA1, DUET_BETA, DUET_OMEGA);
    printf("  gamma1 == B+beta: %s\n", (DUET_GAMMA1==DUET_B+DUET_BETA)?"PASS":"FAIL");
    
    printf("=============================================================\n\n");

    int fail=0;
    printf("  [1/4] Plain Sign/Verify (%d iter)...\n", ITERATIONS);
    if(test_plain())  { printf("  FAIL plain\n");   fail++; } else printf("        PASSED\n");
    printf("  [2/4] Compact Sign/Verify (%d iter)...\n", ITERATIONS);
    if(test_compact()){ printf("  FAIL compact\n"); fail++; } else printf("        PASSED\n");
    printf("  [3/4] rANS Sign/Verify (%d iter)...\n", ITERATIONS);
    if(test_rans())   { printf("  FAIL rans\n");    fail++; } else printf("        PASSED\n");
    printf("  [4/4] Cross-key / wrong-msg / tampered...\n");
    if(test_cross())  { printf("  FAIL cross\n");   fail++; } else printf("        PASSED\n");
    printf("\n");
    if(!fail){
        printf("=============================================================\n");
        printf("  ALL CORRECTNESS TESTS PASSED  (%d iter x 3 variants)\n", ITERATIONS);
        printf("=============================================================\n\n");
    } else {
        printf("=============================================================\n");
        printf("  %d FAILURE(S)\n", fail);
        printf("=============================================================\n\n");
    }
    
    benchmark_cycles();
    return fail;
}
