
#ifndef DUET_SYMMETRIC_H
#define DUET_SYMMETRIC_H

#include "fips202.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef keccak_state xof256_state;
typedef keccak_state stream128_state;
typedef keccak_state stream256_state;

#define STREAM128_BLOCKBYTES  SHAKE128_RATE   /* 168 */
#define STREAM256_BLOCKBYTES  SHAKE256_RATE   /* 136 */

/* SHAKE-128 stream: seed (32 B) || nonce (2 B LE) */
void shake128_stream_init(keccak_state *st,
                          const uint8_t seed[DUET_SEEDBYTES],
                          uint16_t nonce);

/* SHAKE-256 stream: seed (32 B) || nonce (2 B LE) */
void shake256_stream_init(keccak_state *st,
                          const uint8_t seed[DUET_SEEDBYTES],
                          uint16_t nonce);

/* H(in1 || in2) via SHAKE-256 */
void shake256_absorb_twice(keccak_state *st,
                            const uint8_t *in1, size_t in1len,
                            const uint8_t *in2, size_t in2len);

/* H(in1 || in2 || in3) via SHAKE-256 */
void shake256_absorb_thrice(keccak_state *st,
                             const uint8_t *in1, size_t in1len,
                             const uint8_t *in2, size_t in2len,
                             const uint8_t *in3, size_t in3len);

/* Convenience macros  */
#define xof256_absorb_once(S,IN,LEN)  shake256_absorb_once(S,IN,LEN)
#define xof256_squeeze(OUT,LEN,S)     shake256_squeeze(OUT,LEN,S)
#define xof256_absorb_twice(S,I1,L1,I2,L2)  shake256_absorb_twice(S,I1,L1,I2,L2)
#define xof256_absorb_thrice(S,I1,L1,I2,L2,I3,L3) \
        shake256_absorb_thrice(S,I1,L1,I2,L2,I3,L3)

#define stream128_init(S,SEED,NONCE)  shake128_stream_init(S,SEED,NONCE)
#define stream128_squeezeblocks(OUT,NB,S)  shake128_squeezeblocks(OUT,NB,S)
#define stream256_init(S,SEED,NONCE)  shake256_stream_init(S,SEED,NONCE)
#define stream256_squeezeblocks(OUT,NB,S)  shake256_squeezeblocks(OUT,NB,S)

#endif /* !DUET_SYMMETRIC_H */
