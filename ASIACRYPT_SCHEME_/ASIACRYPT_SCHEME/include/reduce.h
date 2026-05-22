
#ifndef DUET_REDUCE_H
#define DUET_REDUCE_H

#include "params.h"
#include <stdint.h>

/*
 * Montgomery constants for q = 64513:
 *   MONT  = 2^32 mod q = 14321
 *   QINV  = q^{-1} mod 2^32 = 940508161
 *   QREC  = floor(2^32 / q) = 66575
 */
#define MONT   14321u
#define QINV   940508161u
#define QREC   66575u

#define montgomery_reduce  DUET_NAMESPACE(montgomery_reduce)
#define barrett_reduce     DUET_NAMESPACE(barrett_reduce)
#define caddq              DUET_NAMESPACE(caddq)
#define freeze             DUET_NAMESPACE(freeze)
#define barrett_reduce2q   DUET_NAMESPACE(barrett_reduce2q)
#define freeze2q           DUET_NAMESPACE(freeze2q)

/* r = a * 2^{-32} mod q */
int32_t montgomery_reduce(int64_t a);
/* narrow a mod q into (-q, q) */
int32_t barrett_reduce(int32_t a);
/* if a < 0, add q */
int32_t caddq(int32_t a);
/* canonical [0, q) */
int32_t freeze(int32_t a);
/* narrow into [0, 2q) approximately */
int32_t barrett_reduce2q(int32_t a);
/* canonical [0, 2q) */
int32_t freeze2q(int32_t a);

#endif /* !DUET_REDUCE_H */
