
#ifndef DUET_ENCODING_H
#define DUET_ENCODING_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

#define encode_signature_rans  DUET_NAMESPACE(encode_signature_rans)
#define decode_signature_rans  DUET_NAMESPACE(decode_signature_rans)

void encode_signature_rans(uint8_t *out, size_t *outlen,
                            const uint8_t *sig, size_t siglen);
int  decode_signature_rans(uint8_t *out, size_t *outlen,
                            const uint8_t *in,  size_t inlen);

#endif
