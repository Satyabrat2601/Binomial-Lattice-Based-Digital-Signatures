
#ifndef DUET_NTT_H
#define DUET_NTT_H

#include "params.h"
#include <stdint.h>

#define ntt           DUET_NAMESPACE(ntt)
#define invntt_tomont DUET_NAMESPACE(invntt_tomont)

/* Forward NTT; output in bit-reversed order */
void ntt(int32_t a[DUET_N]);

/* Inverse NTT multiplied by 2^32 (Montgomery factor) */
void invntt_tomont(int32_t a[DUET_N]);

#endif /* !DUET_NTT_H */
