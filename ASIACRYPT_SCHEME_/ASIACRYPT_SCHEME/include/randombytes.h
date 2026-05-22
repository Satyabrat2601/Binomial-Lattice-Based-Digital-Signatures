
#ifndef DUET_RANDOMBYTES_H
#define DUET_RANDOMBYTES_H

#include <stddef.h>
#include <stdint.h>

/*   Returns 0 on success, -1 on failure. */
int randombytes(uint8_t *buf, size_t outlen);

#endif /* !DUET_RANDOMBYTES_H */
