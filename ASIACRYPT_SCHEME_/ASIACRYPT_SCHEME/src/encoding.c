
#include "encoding.h"
#include "params.h"
#include "rans_byte.h"
#include <string.h>
#include <stdint.h>

#define SCALE_BITS 10u
#define SCALE      (1u << SCALE_BITS)  /* 1024 */

/* ====================================================================
   z1 rANS frequency tables (per-mode, from CBD_{gamma1} convolution).
   ==================================================================== */
#if DUET_CONFIG_MODE == DUET_MODE1
#define Z1_NSYMS 76
#define Z1_ESC   75
static const int16_t z1_vals[75] = {
    -37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,
    -19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,
    5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
    32,33,34,35,36,37
};
static const uint16_t z1_freq[76] = {
    1,1,1,1,1,2,2,2,3,3,4,4,5,6,7,8,9,10,11,12,13,15,16,18,19,21,22,23,25,26,
    27,28,29,30,31,31,31,24,31,31,31,30,29,28,27,26,25,23,22,21,19,18,16,15,13,
    12,11,10,9,8,7,6,5,4,4,3,3,2,2,2,1,1,1,1,1,4
};
static const uint16_t z1_cum[76] = {
    0,1,2,3,4,5,7,9,11,14,17,21,25,30,36,43,51,60,70,81,93,106,121,137,155,174,
    195,217,240,265,291,318,346,375,405,436,467,498,522,553,584,615,645,674,702,
    729,755,780,803,825,846,865,883,899,914,927,939,950,960,969,977,984,990,995,
    999,1003,1006,1009,1011,1013,1015,1016,1017,1018,1019,1020
};

#elif DUET_CONFIG_MODE == DUET_MODE3
#define Z1_NSYMS 88
#define Z1_ESC   87
static const int16_t z1_vals[87] = {
    -43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,
    -24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,
    -3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
    26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43
};
static const uint16_t z1_freq[88] = {
    1,1,1,1,1,1,1,2,2,2,3,3,3,4,4,5,6,6,7,8,9,10,10,11,12,13,14,15,17,18,19,20,
    21,22,22,23,24,25,25,26,26,26,27,25,27,26,26,26,25,25,24,23,22,22,21,20,19,18,
    17,15,14,13,12,11,10,10,9,8,7,6,6,5,4,4,3,3,3,2,2,2,1,1,1,1,1,1,1,5
};
static const uint16_t z1_cum[88] = {
    0,1,2,3,4,5,6,7,9,11,13,16,19,22,26,30,35,41,47,54,62,71,81,91,102,114,127,
    141,156,173,191,210,230,251,273,295,318,342,367,392,418,444,470,497,522,549,
    575,601,627,652,677,701,724,746,768,789,809,828,846,863,878,892,905,917,928,
    938,948,957,965,972,978,984,989,993,997,1000,1003,1006,1008,1010,1012,1013,
    1014,1015,1016,1017,1018,1019
};

#elif DUET_CONFIG_MODE == DUET_MODE5
#define Z1_NSYMS 90
#define Z1_ESC   89
static const int16_t z1_vals[89] = {
    -44,-43,-42,-41,-40,-39,-38,-37,-36,-35,-34,-33,-32,-31,-30,-29,-28,-27,-26,
    -25,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,
    -4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
    25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44
};
static const uint16_t z1_freq[90] = {
    1,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,6,6,7,8,8,9,10,11,12,13,14,15,15,16,17,18,
    19,20,21,22,23,23,24,24,25,25,25,26,23,26,25,25,25,24,24,23,23,22,21,20,19,18,
    17,16,15,15,14,13,12,11,10,9,8,8,7,6,6,5,4,4,3,3,3,2,2,2,2,1,1,1,1,1,1,5
};
static const uint16_t z1_cum[90] = {
    0,1,2,3,4,5,6,8,10,12,14,17,20,23,27,31,36,42,48,55,63,71,80,90,101,113,126,
    140,155,170,186,203,221,240,260,281,303,326,349,373,397,422,447,472,498,521,
    547,572,597,622,646,670,693,716,738,759,779,798,816,833,849,864,879,893,906,
    918,929,939,948,956,964,971,977,983,988,992,996,999,1002,1005,1007,1009,1011,
    1013,1014,1015,1016,1017,1018,1019
};
#endif

/* Decode LUT */
static uint8_t z1_lut[SCALE];
static int     z1_lut_built = 0;

static void build_z1_lut(void) {
    if (z1_lut_built) return;
    for (int s = 0; s < Z1_NSYMS; ++s) {
        uint32_t end = (s+1 < Z1_NSYMS) ? z1_cum[s+1] : SCALE;
        for (uint32_t j = z1_cum[s]; j < end; ++j) z1_lut[j] = (uint8_t)s;
    }
    z1_lut_built = 1;
}

static inline int coeff_to_sym(int32_t v) {
    int lo=0, hi=Z1_ESC-1;
    while(lo<=hi){ int m=(lo+hi)>>1; if(z1_vals[m]<v) lo=m+1; else if(z1_vals[m]>v) hi=m-1; else return m; }
    return Z1_ESC;
}

/* Frequencies must sum to SCALE=1024. Empirical measurements: */
#if DUET_CONFIG_MODE == DUET_MODE1
/* p0~0.967 p1~0.018 p2~0.015: freq0=990 freq1=18 freq2=16 */
#define H_FREQ0  990u
#define H_FREQ1   18u
#define H_FREQ2   16u
#elif DUET_CONFIG_MODE == DUET_MODE3
/* slightly different due to k=3 */
#define H_FREQ0  988u
#define H_FREQ1   19u
#define H_FREQ2   17u
#elif DUET_CONFIG_MODE == DUET_MODE5
/* k=4 */
#define H_FREQ0  987u
#define H_FREQ1   19u
#define H_FREQ2   18u
#endif
/* Cumulative: sym0: [0,H_FREQ0), sym1: [H_FREQ0, H_FREQ0+H_FREQ1), sym2: [...,1024) */
#define H_CUM0  0u
#define H_CUM1  H_FREQ0
#define H_CUM2  (H_FREQ0 + H_FREQ1)


void encode_signature_rans(uint8_t *out, size_t *outlen,
                           const uint8_t *sig, size_t siglen)
{
    if (!out || !outlen || !sig || siglen < (size_t)DUET_CRYPTO_BYTES) {
        if (outlen) *outlen = 0; return;
    }
    build_z1_lut();

    const unsigned z1_n     = (DUET_ELL + 1) * DUET_N;
    const unsigned z1_raw_b = (DUET_ELL + 1) * DUET_POLYZ_PACKEDBYTES;
    const unsigned h_raw_b  = DUET_K * DUET_POLYH_PACKEDBYTES;

    /* --- Unpack z1 from 10-bit packed format --- */
    int32_t z1v[z1_n];
    {
        const uint8_t *zp = sig;
        uint64_t acc=0; int bits=0; unsigned pos=0;
        const uint32_t mask=(1u<<DUET_B_BITS)-1u;
        for (unsigned i=0; i<z1_n; ++i) {
            while(bits<(int)DUET_B_BITS){acc|=((uint64_t)zp[pos++]<<bits);bits+=8;}
            z1v[i]=(int32_t)(acc&mask)-(int32_t)DUET_B;
            acc>>=DUET_B_BITS; bits-=(int)DUET_B_BITS;
        }
    }

    const uint8_t *h_raw = sig + z1_raw_b;

    /* --- Extract c_seed --- */
    const uint8_t *c_seed = sig + z1_raw_b + h_raw_b;

    /* --- Build ESC tail bitstream for z1 --- */
    /* Count tail values and build tail bitstream */
    unsigned tail_cnt = 0;
    for (unsigned i=0; i<z1_n; ++i) if(coeff_to_sym(z1v[i])==Z1_ESC) tail_cnt++;

    uint8_t tail_bits[(z1_n*10/8)+4];
    unsigned tail_bit_off=0;
    memset(tail_bits,0,sizeof(tail_bits));
    for (unsigned i=0; i<z1_n; ++i) {
        if(coeff_to_sym(z1v[i])!=Z1_ESC) continue;
        uint32_t raw=(uint32_t)(z1v[i]+(int32_t)DUET_B);
        for(unsigned b=0;b<10u;++b){
            if(raw&(1u<<b)) tail_bits[tail_bit_off>>3]|=(uint8_t)(1u<<(tail_bit_off&7));
            ++tail_bit_off;
        }
    }
    unsigned tail_bytes=(tail_bit_off+7)/8;

    /* --- Joint rANS stream encoding (backwards order): h, then z1 --- */
    /* Buffer sized for z1 + h bits + overhead */
    uint8_t rans_buf[z1_n*2 + DUET_K*DUET_N/8 + 64];
    uint8_t *rans_end = rans_buf + sizeof(rans_buf);
    uint8_t *rans_ptr = rans_end;
    RansState R; RansEncInit(&R);

    /* Encode h (ternary) backwards: n*k symbols.
     * sym=0: h==0 (no correction), sym=1: h==1 (+1), sym=2: h==H_RANGE-1 (-1 mod H_RANGE) */
    for (int i=(int)(DUET_K*DUET_N)-1; i>=0; --i) {
        uint8_t hv = h_raw[i];
        int sym = (hv == 0) ? 0 : (hv == 1) ? 1 : 2;
        uint32_t cums[3]  = {H_CUM0, H_CUM1, H_CUM2};
        uint32_t freqs[3] = {H_FREQ0, H_FREQ1, H_FREQ2};
        RansEncPut(&R, &rans_ptr, cums[sym], freqs[sym], SCALE_BITS);
    }

    /* Encode z1 symbols backwards */
    for (int i=(int)z1_n-1; i>=0; --i) {
        int s=coeff_to_sym(z1v[i]);
        RansEncPut(&R, &rans_ptr, z1_cum[s], z1_freq[s], SCALE_BITS);
    }
    RansEncFlush(&R, &rans_ptr);

    size_t rans_bytes = (size_t)(rans_end - rans_ptr);
    size_t z_len = rans_bytes + tail_bytes;

    /* Output: [2B z_len] [rANS stream] [tail bits] [c_seed 32B] */
    size_t op=0;
    out[op++]=(uint8_t)(z_len);
    out[op++]=(uint8_t)(z_len>>8);
    memcpy(out+op, rans_ptr, rans_bytes); op+=rans_bytes;
    memcpy(out+op, tail_bits, tail_bytes); op+=tail_bytes;
    /* c_seed: store verbatim (32B) - close to 25B entropy target */
    memcpy(out+op, c_seed, DUET_POLYC_PACKEDBYTES); op+=DUET_POLYC_PACKEDBYTES;
    *outlen=op;
}

/* ====================================================================
   decode_signature_rans
   ==================================================================== */
int decode_signature_rans(uint8_t *out, size_t *outlen,
                          const uint8_t *in, size_t inlen)
{
    if (!out||!outlen||!in) return -1;
    build_z1_lut();
    if (inlen < 2u+4u+DUET_POLYC_PACKEDBYTES) return -1;

    size_t z_len=(size_t)in[0]|((size_t)in[1]<<8);
    if (inlen < 2u+z_len+DUET_POLYC_PACKEDBYTES) return -1;

    const uint8_t *z_block  = in+2;
    const uint8_t *c_seed_p = in+2+z_len;
    const unsigned z1_n=(DUET_ELL+1)*DUET_N;

    /* Decode rANS stream: z1 symbols first, then h bits */
    uint8_t *rp=(uint8_t*)z_block;
    const uint8_t *rp_end=z_block+z_len;
    RansState Rs; if(RansDecInit(&Rs,&rp)) return -1;

    int32_t z1v[z1_n]; int syms[z1_n]; unsigned tail_cnt=0;
    for (unsigned i=0; i<z1_n; ++i) {
        uint32_t cf=RansDecGet(&Rs,SCALE_BITS);
        int s=(int)z1_lut[cf]; syms[i]=s;
        if(s==Z1_ESC) ++tail_cnt;
        RansDecAdvance(&Rs,&rp,rp_end,z1_cum[s],z1_freq[s],SCALE_BITS);
    }

    /* Decode h (ternary): restore actual correction values.
     * sym=0 -> h=0, sym=1 -> h=1, sym=2 -> h=H_RANGE-1 (-1 correction) */
    uint8_t h_vals[DUET_K*DUET_N];
    for (unsigned i=0; i<DUET_K*DUET_N; ++i) {
        uint32_t cf=RansDecGet(&Rs,SCALE_BITS);
        int sym;
        uint32_t cum, freq;
        if      (cf < H_CUM1) { sym=0; cum=H_CUM0; freq=H_FREQ0; }
        else if (cf < H_CUM2) { sym=1; cum=H_CUM1; freq=H_FREQ1; }
        else                  { sym=2; cum=H_CUM2; freq=H_FREQ2; }
        /* Restore actual h value: sym=2 maps to H_RANGE-1 */
        h_vals[i] = (sym==0) ? 0 : (sym==1) ? 1 : (uint8_t)(DUET_H_RANGE - 1);
        RansDecAdvance(&Rs,&rp,rp_end,cum,freq,SCALE_BITS);
    }

    /* Reconstruct z1 values from symbols + tail */
    const uint8_t *tail_p=rp; unsigned tboff=0;
    for (unsigned i=0; i<z1_n; ++i) {
        int s=syms[i];
        if(s<Z1_ESC) { z1v[i]=(int32_t)z1_vals[s]; continue; }
        uint32_t raw=0;
        for(unsigned b=0;b<10u;++b){
            if(tail_p[tboff>>3]&(1u<<(tboff&7))) raw|=(1u<<b); ++tboff;
        }
        z1v[i]=(int32_t)raw-(int32_t)DUET_B;
    }

    /* Pack z1 as 10-bit */
    const unsigned z1_raw_b=(DUET_ELL+1)*DUET_POLYZ_PACKEDBYTES;
    uint8_t *op=out;
    {
        uint64_t acc=0; int bits=0;
        for(unsigned i=0;i<z1_n;++i){
            int32_t v=z1v[i];
            if(v<-(int32_t)DUET_B) v=-(int32_t)DUET_B;
            if(v> (int32_t)DUET_B) v= (int32_t)DUET_B;
            uint32_t u=(uint32_t)(v+(int32_t)DUET_B);
            acc|=((uint64_t)u<<bits); bits+=(int)DUET_B_BITS;
            while(bits>=8){*op++=(uint8_t)acc;acc>>=8;bits-=8;}
        }
        if(bits>0)*op++=(uint8_t)acc;
        while((size_t)(op-out)<z1_raw_b)*op++=0;
        op=out+z1_raw_b;
    }

    /* Restore h with actual correction values (0, 1, or H_RANGE-1) */
    for(unsigned i=0;i<DUET_K;++i){
        for(unsigned r=0;r<DUET_N;++r)
            op[r]=h_vals[i*DUET_N+r];
        op+=DUET_POLYH_PACKEDBYTES;
    }

    /* c_seed */
    memcpy(op, c_seed_p, DUET_POLYC_PACKEDBYTES); op+=DUET_POLYC_PACKEDBYTES;
    *outlen=(size_t)(op-out);
    return 0;
}
