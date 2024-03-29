#ifndef _HWBSPECK_H_
#define _HWBSPECK_H_
#include "WBMatrix/WBMatrix.h"
#include "math.h"
#include "string.h"

#define u8 uint8_t
#define u16 uint16_t
#define u32 uint32_t
#define u64 uint64_t

#define SPECK_ROUNDS 22
#define SPECK_KEY_LEN 4

#define pt 16
#define apt 256
#define tb32 32
#define tb64 64
#define tb128 128
static const u16 key16[4] = {0x1918, 0x1110, 0x0908, 0x0100};
static const u32 key32[4] = {0x1b1a1918, 0x13121110, 0x0b0a0908, 0x03020100};
static const u64 key64[4] = {0x1f1e1d1c1b1a1918, 0x1716151413121110, 0x0f0e0d0c0b0a0908, 0x0706050403020100};
static const int ktb = 0x0;
static const int kte = 0xf;
#define R32(x, y, k) (x = ROR32(x, 7), x += y, x ^= k, y = ROL32(y, 2), y ^= x)
#define R64(x, y, k) (x = ROR64(x, 8), x += y, x ^= k, y = ROL64(y, 3), y ^= x)
#define R128(x, y, k) (x = ROR128(x, 8), x += y, x ^= k, y = ROL128(y, 3), y ^= x)

static const double EPS = 1e-6;
void Key_expand32(u16 key[4], u16 roundkey[SPECK_ROUNDS]);
void Key_expand64(u32 key[4], u32 roundkey[SPECK_ROUNDS]);
void Key_expand128(u64 key[4], u64 roundkey[SPECK_ROUNDS]);

void DCA32();
void sub_DCA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff, int q, u16 *reduced_key, int *reduced_key_count);
void DCA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff);

void DCA64();
void sub_DCA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff, int q, u32 *reduced_key, int *reduced_key_count);
void DCA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff);

void DCA128();
void sub_DCA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff, int q, u64 *reduced_key, int *reduced_key_count);
void DCA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff);

u16 ROL32(u16 x, int r);
u16 ROR32(u16 x, int r);
u32 ROL64(u32 x, int r);
u32 ROR64(u32 x, int r);
u64 ROL128(u64 x, int r);
u64 ROR128(u64 x, int r);
#endif