#include "speed_print.h"
#include "cpucycles.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int cmp_u64(const void *a, const void *b) {
    if (*(uint64_t*)a < *(uint64_t*)b) return -1;
    if (*(uint64_t*)a > *(uint64_t*)b) return  1;
    return 0;
}
static uint64_t median(uint64_t *l, size_t n) {
    qsort(l, n, sizeof(uint64_t), cmp_u64);
    return (n & 1) ? l[n/2] : (l[n/2-1] + l[n/2]) / 2;
}
static uint64_t average(uint64_t *t, size_t n) {
    uint64_t s = 0;
    for (size_t i = 0; i < n; ++i) s += t[i];
    return s / n;
}

void print_results(const char *s, uint64_t *t, size_t tlen) {
    static uint64_t overhead = UINT64_MAX;
    if (overhead == UINT64_MAX) overhead = cpucycles_overhead();
    if (tlen < 2) { fprintf(stderr, "Need >= 2 samples\n"); return; }
    --tlen;
    for (size_t i = 0; i < tlen; ++i)
        t[i] = t[i+1] - t[i] - overhead;
    printf("%-44s  median: %8llu  avg: %8llu  cycles\n",
           s,
           (unsigned long long)median(t, tlen),
           (unsigned long long)average(t, tlen));
}
