
#include "randombytes.h"
#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <unistd.h>

int randombytes(uint8_t *buf, size_t outlen) {
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) return -1;
    while (outlen > 0) {
        ssize_t r = read(fd, buf, outlen);
        if (r < 0) {
            if (errno == EINTR) continue;
            close(fd);
            return -1;
        }
        buf    += (size_t)r;
        outlen -= (size_t)r;
    }
    close(fd);
    return 0;
}
