
#ifndef DUET_CPUCYCLES_H
#define DUET_CPUCYCLES_H

#include <stdint.h>

#if defined(__x86_64__) || defined(__i386__)
static inline uint64_t cpucycles(void) {
    uint64_t result;
    __asm__ volatile("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
                     : "=a"(result) : : "%rdx");
    return result;
}
#elif defined(__aarch64__)
static inline uint64_t cpucycles(void) {
    uint64_t result;
    __asm__ volatile("mrs %0, cntvct_el0" : "=r"(result));
    return result;
}
#else
#include <time.h>
static inline uint64_t cpucycles(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}
#endif

uint64_t cpucycles_overhead(void);

#endif /* !DUET_CPUCYCLES_H */
