
#include "reduce.h"
#include <stdint.h>

int32_t montgomery_reduce(int64_t a) {
    int32_t t = (int32_t)((uint32_t)a * QINV);
    t = (int32_t)((a - (int64_t)t * DUET_Q) >> 32);
    return t;
}

int32_t barrett_reduce(int32_t a) {
    int32_t t = (int32_t)(((int64_t)a * QREC) >> 32);
    return a - t * (int32_t)DUET_Q;
}

int32_t caddq(int32_t a) {
    a += (a >> 31) & (int32_t)DUET_Q;
    return a;
}

int32_t freeze(int32_t a) {
    a = barrett_reduce(a);
    a = caddq(a);
    return a;
}

int32_t freeze2q(int32_t a) {
    int32_t dq = (int32_t)DUET_DQ;
    a = (int32_t)(((int64_t)a % dq + dq) % dq);
    return a;
}

int32_t barrett_reduce2q(int32_t a) {
    return freeze2q(a);
}
