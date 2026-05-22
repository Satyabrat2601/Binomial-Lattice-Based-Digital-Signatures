
#ifndef DUET_PACKING_H
#define DUET_PACKING_H
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <stdint.h>

#define pack_pk    DUET_NAMESPACE(pack_pk)
#define unpack_pk  DUET_NAMESPACE(unpack_pk)
void pack_pk(uint8_t pk[DUET_CRYPTO_PUBLICKEYBYTES],
             const uint8_t rho[DUET_SEEDBYTES], const polyveck *b);
void unpack_pk(uint8_t rho[DUET_SEEDBYTES], polyveck *b,
               const uint8_t pk[DUET_CRYPTO_PUBLICKEYBYTES]);

#define pack_sk    DUET_NAMESPACE(pack_sk)
#define unpack_sk  DUET_NAMESPACE(unpack_sk)
void pack_sk(uint8_t sk[DUET_CRYPTO_SECRETKEYBYTES],
             const uint8_t rho[DUET_SEEDBYTES], const uint8_t key[DUET_SEEDBYTES],
             const uint8_t tr[DUET_SEEDBYTES],
             const polyvecell *s1, const polyveck *s2);
void unpack_sk(uint8_t rho[DUET_SEEDBYTES], uint8_t key[DUET_SEEDBYTES],
               uint8_t tr[DUET_SEEDBYTES], polyvecell *s1, polyveck *s2,
               const uint8_t sk[DUET_CRYPTO_SECRETKEYBYTES]);

#define pack_sig    DUET_NAMESPACE(pack_sig)
#define unpack_sig  DUET_NAMESPACE(unpack_sig)
void pack_sig(uint8_t sig[DUET_CRYPTO_BYTES],
              const polyvecL *z1, const polyveck *h,
              const uint8_t c_seed[DUET_SEEDBYTES]);
int  unpack_sig(polyvecL *z1, polyveck *h, uint8_t c_seed[DUET_SEEDBYTES],
                const uint8_t sig[DUET_CRYPTO_BYTES]);

#endif /* DUET_PACKING_H */
