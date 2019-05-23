/*
  Deflicker plugin for Avisynth 2.6 -  mean intensity stabilizer
  Version 0.4 August 16,2004
    (c) 2004, A.G.Balakhnin aka Fizick
  bag@hotmail.ru

  Code update, AVS2.6, x64, SSE2, AVX2 2019 by pinterf

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


*/
//following  includes needed
#include "windows.h"
#include "avisynth.h"

#include "stdio.h"  // for debug
#include "math.h"  // for sqrt, fabs

#include "info.h" // for info on frame
#include "stdint.h"
#include "immintrin.h"


template<bool isYUY2>
static __forceinline void SumLine_avx2(const BYTE* srcp, int width, int* mean, int* varline)
{
  const int wmod32 = width / 32 * 32;
  auto zero = _mm256_setzero_si256();
  int sum = 0;
  int sqrsum = 0;
  auto sum_simd = _mm256_setzero_si256();
  auto sqrsum_simd = _mm256_setzero_si256();
  for (int x = 0; x < wmod32; x += 32)
  {
    auto src = _mm256_loadu_si256((const __m256i*)(srcp + x));
    if constexpr (isYUY2)
      src = _mm256_and_si256(src, _mm256_set1_epi16(0x00FF)); // luma mask
    auto horiz_sum = _mm256_sad_epu8(src, zero);
    sum_simd = _mm256_add_epi32(sum_simd, horiz_sum);

    auto src_lo = _mm256_unpacklo_epi8(src, zero);
    auto mul_lo = _mm256_mullo_epi16(src_lo, src_lo);
    auto src_hi = _mm256_unpackhi_epi8(src, zero);
    auto mul_hi = _mm256_mullo_epi16(src_hi, src_hi);

    // A = a*a
    // Ahi Alo Ahi Alo Ahi Alo Ahi Alo Ahi Alo Ahi Alo Ahi Alo Ahi Alo
    // Alo   0 Alo   0 Alo   0 Alo   0 Alo   0 Alo   0 Alo   0 Alo   0
    // Alo+0+Alo+0+Alo+0+Alo+0         Alo+0+Alo+0+Alo+0+Alo+0
    auto horiz_sum16_lo_lo8 = _mm256_sad_epu8(_mm256_slli_epi16(mul_lo, 8), zero); // sum LSB
    sqrsum_simd = _mm256_add_epi32(sqrsum_simd, horiz_sum16_lo_lo8);
    //   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0
    // 0+Ahi+0+Ahi+0+Ahi+0+Ahi           0+Ahi+0+Ahi+0+Ahi+0+Ahi
    auto horiz_sum16_lo_hi8 = _mm256_sad_epu8(_mm256_srli_epi16(mul_lo, 8), zero); // sum MSB, later to be multiplied by 256
    sqrsum_simd = _mm256_add_epi32(sqrsum_simd, _mm256_slli_epi32(horiz_sum16_lo_hi8, 8));

    auto horiz_sum16_hi_lo8 = _mm256_sad_epu8(_mm256_slli_epi16(mul_hi, 8), zero); // sum LSB
    sqrsum_simd = _mm256_add_epi32(sqrsum_simd, horiz_sum16_hi_lo8);
    //   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0 Ahi   0
    // 0+Ahi+0+Ahi+0+Ahi+0+Ahi           0+Ahi+0+Ahi+0+Ahi+0+Ahi
    auto horiz_sum16_hi_hi8 = _mm256_sad_epu8(_mm256_srli_epi16(mul_hi, 8), zero); // sum MSB, later to be multiplied by 256
    sqrsum_simd = _mm256_add_epi32(sqrsum_simd, _mm256_slli_epi32(horiz_sum16_hi_hi8, 8));

  }

  // sum here: two 32 bit partial result: sum1 0 sum2 0
  __m128i sum_simd_128 = _mm_add_epi32(_mm256_castsi256_si128(sum_simd), _mm256_extractf128_si256(sum_simd, 1));
  auto sum_hi = _mm_unpackhi_epi64(sum_simd_128, _mm256_castsi256_si128(zero));
  sum_simd_128 = _mm_add_epi32(sum_simd_128, sum_hi);
  int rowsum = _mm_cvtsi128_si32(sum_simd_128);
  *mean += rowsum;

  // sum here: two 32 bit partial result: sum1 0 sum2 0
  __m128i sqrsum_simd_128 = _mm_add_epi32(_mm256_castsi256_si128(sqrsum_simd), _mm256_extractf128_si256(sqrsum_simd, 1));
  auto sqrsum_hi = _mm_unpackhi_epi64(sqrsum_simd_128, _mm256_castsi256_si128(zero));
  sqrsum_simd_128 = _mm_add_epi32(sqrsum_simd_128, sqrsum_hi);
  int rowsqrsum = _mm_cvtsi128_si32(sqrsum_simd_128);
  *varline += rowsqrsum;

  int increment;
  if constexpr (isYUY2)
    increment = 2;
  else
    increment = 1;


  for (int x = wmod32; x < width; x += increment)
  {   // don't use step=2 for non YUY2
    int cur = srcp[x];
    *mean += cur;
    *varline += cur * cur; // simply quadrat, sub mean later
  }

}

static __forceinline void CorrectYUV_avx2(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 

  const int wmod32 = src_width / 32 * 32;
  auto zero = _mm256_setzero_si256();
  const __m256i mult256w_simd = _mm256_set1_epi16(mult256w);
  const __m256i addw_simd = _mm256_set1_epi16(addw);
  const __m256i lmin_simd = _mm256_set1_epi8(lmin);
  const __m256i lmax_simd = _mm256_set1_epi8(lmax);
  for (int x = 0; x < wmod32; x += 32)
  {
    auto src = _mm256_loadu_si256((const __m256i*)(srcp + x));
    auto src_lo_mul256 = _mm256_unpacklo_epi8(zero, src);
    auto src_hi_mul256 = _mm256_unpackhi_epi8(zero, src);

    auto mul_lo = _mm256_mulhi_epu16(src_lo_mul256, mult256w_simd);
    auto mul_hi = _mm256_mulhi_epu16(src_hi_mul256, mult256w_simd);

    auto tmp_lo = _mm256_adds_epi16(mul_lo, addw_simd);
    auto tmp_hi = _mm256_adds_epi16(mul_hi, addw_simd);

    auto res = _mm256_min_epu8(_mm256_max_epu8(_mm256_packus_epi16(tmp_lo, tmp_hi), lmin_simd), lmax_simd);
    _mm256_storeu_si256((__m256i*)(dstp + x), res);
  }

  for (int x = wmod32; x < src_width; x++)
  {
    int cur = ((srcp[x] * mult256w) >> 8) + addw;
    dstp[x] = min(max(cur, lmin), lmax);
  }

}

static __forceinline void CorrectYUY2_avx2(BYTE* dstp, const BYTE* srcp, int src_width, short mult256w, short addw, short lmin, short lmax)
{ // correct luma 
  const int wmod32 = src_width / 32 * 32;
  auto zero = _mm256_setzero_si256();
  const __m256i mult256w_simd = _mm256_set1_epi16(mult256w);
  const __m256i addw_simd = _mm256_set1_epi16(addw);
  const __m256i lmin_simd = _mm256_set1_epi8(lmin);
  const __m256i lmax_simd = _mm256_set1_epi8(lmax);
  const __m256i uvmask = _mm256_set1_epi16(0xFF00);
  for (int x = 0; x < wmod32; x += 32)
  {
    auto src = _mm256_loadu_si256((const __m256i*)(srcp + x)); // YUYV
    auto uv = _mm256_and_si256(src, uvmask);
    auto src_lo_mul256 = _mm256_unpacklo_epi8(zero, src);
    auto src_hi_mul256 = _mm256_unpackhi_epi8(zero, src);

    auto mul_lo = _mm256_mulhi_epu16(src_lo_mul256, mult256w_simd);
    auto mul_hi = _mm256_mulhi_epu16(src_hi_mul256, mult256w_simd);

    auto tmp_lo = _mm256_adds_epi16(mul_lo, addw_simd);
    auto tmp_hi = _mm256_adds_epi16(mul_hi, addw_simd);

    auto res = _mm256_min_epu8(_mm256_max_epu8(_mm256_packus_epi16(tmp_lo, tmp_hi), lmin_simd), lmax_simd);
    res = _mm256_andnot_si256(uvmask, res);
    res = _mm256_or_si256(res, uv);
    _mm256_storeu_si256((__m256i*)(dstp + x), res);
  }

  for (int x = wmod32; x < src_width; x += 2)
  {
    int cur = ((srcp[x] * mult256w) >> 8) + addw;
    dstp[x] = min(max(cur, lmin), lmax); //Y
    dstp[x + 1] = srcp[x + 1]; // U,V
  }

}

void SumFrame_YUY2_avx2(const BYTE* srcp, int src_pitch, int width, int height, int borderw, int border, int64_t* mean64, int64_t* var64)
{
  for (int h = border; h < height - border; h++)
  {   // Loop from top line to bottom line 
    int varline = 0; // line luma variation 
    int meanline = 0;
    SumLine_avx2<true>(srcp, width - borderw * 2, &meanline, &varline);
    *var64 += varline;
    *mean64 += meanline;
    srcp += src_pitch;
  }
  _mm256_zeroupper();
}

void SumFrame_avx2(const BYTE* srcp, int src_pitch, int width, int height, int borderw, int border, int64_t* mean64, int64_t* var64)
{
  for (int h = border; h < height - border; h++)
  {   // Loop from top line to bottom line 
    int varline = 0; // line luma variation 
    int meanline = 0;
    SumLine_avx2<false>(srcp, width - borderw * 2, &meanline, &varline);
    *var64 += varline;
    *mean64 += meanline;
    srcp += src_pitch;
  }
  _mm256_zeroupper();
}

void CorrectFrame_YUY2_avx2(BYTE* dstp, int dst_pitch, const BYTE* srcp, int src_pitch, int src_height, int src_width, short mult256w, short addw, short lmin, short lmax)
{
  for (int h = 0; h < src_height; h++)
  {       // Loop from top line to bottom line 
    CorrectYUY2_avx2(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
    srcp += src_pitch;
    dstp += dst_pitch;
  }
  _mm256_zeroupper();
}

void CorrectFrame_avx2(BYTE* dstp, int dst_pitch, const BYTE* srcp, int src_pitch, int src_height, int src_width, short mult256w, short addw, short lmin, short lmax)
{
  for (int h = 0; h < src_height; h++)
  {       // Loop from top line to bottom line 
    CorrectYUV_avx2(dstp, srcp, src_width, mult256w, addw, lmin, lmax);
    srcp += src_pitch;
    dstp += dst_pitch;
  }
  _mm256_zeroupper();
}