
#ifndef RANS_BYTE_HEADER
#define RANS_BYTE_HEADER

#include <stdint.h>

#ifdef assert
#define RansAssert assert
#else
#define RansAssert(x)
#endif

#define RANS_BYTE_L (1u << 23) // lower bound of our normalization interval

typedef uint32_t RansState;

// Initialize a rANS encoder.
static inline void RansEncInit(RansState *r) { *r = RANS_BYTE_L; }

// Renormalize the encoder. Internal function.
static inline RansState RansEncRenorm(RansState x, uint8_t **pptr,
                                      uint32_t freq, uint32_t scale_bits) {
  uint32_t x_max =
      ((RANS_BYTE_L >> scale_bits) << 8) * freq; // this turns into a shift.
  if (x >= x_max) {
    uint8_t *ptr = *pptr;
    do {
      *--ptr = (uint8_t)(x & 0xff);
      x >>= 8;
    } while (x >= x_max);
    *pptr = ptr;
  }
  return x;
}

static inline void RansEncPut(RansState *r, uint8_t **pptr, uint32_t start,
                              uint32_t freq, uint32_t scale_bits) {
  // renormalize
  RansState x = RansEncRenorm(*r, pptr, freq, scale_bits);

  // x = C(s,x)
  *r = ((x / freq) << scale_bits) + (x % freq) + start;
}

// Flushes the rANS encoder.
static inline void RansEncFlush(RansState *r, uint8_t **pptr) {
  uint32_t x = *r;
  uint8_t *ptr = *pptr;

  ptr -= 4;
  ptr[0] = (uint8_t)(x >> 0);
  ptr[1] = (uint8_t)(x >> 8);
  ptr[2] = (uint8_t)(x >> 16);
  ptr[3] = (uint8_t)(x >> 24);

  *pptr = ptr;
}

// Initializes a rANS decoder.
// Unlike the encoder, the decoder works forwards as you'd expect.
static inline int RansDecInit(RansState *r, uint8_t **pptr) {
  uint32_t x;
  uint8_t *ptr = *pptr;

  x = (uint32_t)ptr[0] << 0;
  x |= (uint32_t)ptr[1] << 8;
  x |= (uint32_t)ptr[2] << 16;
  x |= (uint32_t)ptr[3] << 24;
  if (x < RANS_BYTE_L || (RANS_BYTE_L << 8) <= x)
    return 1; // initial state out of range

  ptr += 4;
  *pptr = ptr;
  *r = x;
  return 0;
}

// Returns the current cumulative frequency (map it to a symbol yourself!)
static inline uint32_t RansDecGet(RansState *r, uint32_t scale_bits) {
  return *r & ((1u << scale_bits) - 1);
}

static inline void RansDecAdvance(RansState *r, uint8_t **pptr,
                                  const uint8_t *end, uint32_t start,
                                  uint32_t freq, uint32_t scale_bits) {
  uint32_t mask = (1u << scale_bits) - 1;

  // s, x = D(x)
  uint32_t x = *r;
  x = freq * (x >> scale_bits) + (x & mask) - start;

  // renormalize
  if (x < RANS_BYTE_L && *pptr < end) {
    uint8_t *ptr = *pptr;
    do
      x = (x << 8) | *ptr++;
    while (x < RANS_BYTE_L && ptr < end);
    *pptr = ptr;
  }

  *r = x;
}

typedef struct {
  uint32_t x_max;     // (Exclusive) upper bound of pre-normalization interval
  uint32_t rcp_freq;  // Fixed-point reciprocal frequency
  uint32_t bias;      // Bias
  uint16_t cmpl_freq; // Complement of frequency: (1 << scale_bits) - freq
  uint16_t rcp_shift; // Reciprocal shift
} RansEncSymbol;

// Decoder symbols are straightforward.
typedef struct {
  uint16_t start; // Start of range.
  uint16_t freq;  // Symbol frequency.
} RansDecSymbol;

// Initializes an encoder symbol to start "start" and frequency "freq"
static inline void RansEncSymbolInit(RansEncSymbol *s, uint32_t start,
                                     uint32_t freq, uint32_t scale_bits) {
  RansAssert(scale_bits <= 16);
  RansAssert(start <= (1u << scale_bits));
  RansAssert(freq <= (1u << scale_bits) - start);


  s->x_max = ((RANS_BYTE_L >> scale_bits) << 8) * freq;
  s->cmpl_freq = (uint16_t)((1 << scale_bits) - freq);
  if (freq < 2) {
    s->rcp_freq = ~0u;
    s->rcp_shift = 0;
    s->bias = start + (1 << scale_bits) - 1;
  } else {
    // Alverson, "Integer Division using reciprocals"
    // shift=ceil(log2(freq))
    uint32_t shift = 0;
    while (freq > (1u << shift))
      shift++;

    s->rcp_freq = (uint32_t)(((1ull << (shift + 31)) + freq - 1) / freq);
    s->rcp_shift = shift - 1;

    // With these values, 'q' is the correct quotient, so we
    // have bias=start.
    s->bias = start;
  }
}

// Initialize a decoder symbol to start "start" and frequency "freq"
static inline void RansDecSymbolInit(RansDecSymbol *s, uint32_t start,
                                     uint32_t freq) {
  RansAssert(start <= (1 << 16));
  RansAssert(freq <= (1 << 16) - start);
  s->start = (uint16_t)start;
  s->freq = (uint16_t)freq;
}

static inline void RansEncPutSymbol(RansState *r, uint8_t **pptr,
                                    RansEncSymbol const *sym) {
  RansAssert(sym->x_max != 0); // can't encode symbol with freq=0

  // renormalize
  uint32_t x = *r;
  uint32_t x_max = sym->x_max;
  if (x >= x_max) {
    uint8_t *ptr = *pptr;
    do {
      *--ptr = (uint8_t)(x & 0xff);
      x >>= 8;
    } while (x >= x_max);
    *pptr = ptr;
  }

  uint32_t q =
      (uint32_t)(((uint64_t)x * sym->rcp_freq) >> 32) >> sym->rcp_shift;
  *r = x + sym->bias + q * sym->cmpl_freq;
}

static inline void RansDecAdvanceSymbol(RansState *r, uint8_t **pptr,
                                        const uint8_t *end,
                                        RansDecSymbol const *sym,
                                        uint32_t scale_bits) {
  RansDecAdvance(r, pptr, end, sym->start, sym->freq, scale_bits);
}

static inline void RansDecAdvanceStep(RansState *r, uint32_t start,
                                      uint32_t freq, uint32_t scale_bits) {
  uint32_t mask = (1u << scale_bits) - 1;

  // s, x = D(x)
  uint32_t x = *r;
  *r = freq * (x >> scale_bits) + (x & mask) - start;
}

// Equivalent to RansDecAdvanceStep that takes a symbol.
static inline void RansDecAdvanceSymbolStep(RansState *r,
                                            RansDecSymbol const *sym,
                                            uint32_t scale_bits) {
  RansDecAdvanceStep(r, sym->start, sym->freq, scale_bits);
}

// Renormalize.
static inline void RansDecRenorm(RansState *r, uint8_t **pptr) {
  // renormalize
  uint32_t x = *r;
  if (x < RANS_BYTE_L) {
    uint8_t *ptr = *pptr;
    do
      x = (x << 8) | *ptr++;
    while (x < RANS_BYTE_L);
    *pptr = ptr;
  }

  *r = x;
}

// Verify final state
static inline int RansDecVerify(const RansState *const r) {
  if (*r != RANS_BYTE_L) {
    return 1; // the final state is inconsistent with the initial state
  }
  return 0;
}

#endif // RANS_BYTE_HEADER
